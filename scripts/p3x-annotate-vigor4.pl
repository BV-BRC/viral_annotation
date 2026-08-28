#
# Use vigor4 to annotate the given genome.
#


use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use gjogenbank;
use gjoseqlib;
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use File::Copy;
use IPC::Run qw(run);
use File::SearchPath qw(searchpath);
use Bio::BVBRC::ViralAnnotation::Config qw(vigor_reference_db_directory);
use Bio::BVBRC::ViralAnnotation::VigorTaxonMap;
use P3DataAPI;

use Cwd;

my($opt, $usage) = describe_options("%c %o",
				    ["reference=s" => "Vigor4 reference name"],
				    ["taxon=i" => "Taxon identifier"],
				    ["remove-existing" => "Remove existing CDS and mat_peptide features if vigor4 run is successful"],	
				    ["input|i=s" => "Input file"],
				    ["output|o=s" => "Output file"],
				    ["threads|j=i" => "Limit vigor to this many threads", { default => 1 }],
				    ["debug|d" => "Enable debugging"],
				    ["help|h" => "Show this help message"]);
print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;

chomp(my $hostname = `hostname`);

my $tempdir = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1));

print STDERR "Tempdir=$tempdir\n" if $opt->debug;

my $here = getcwd;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";

#
# Determine our reference database.
#
# If --reference passed, use that.
#
# Otherwise get the taxon id from --taxon parameter or from the GTO.
#

#
# Metadata accumulated for the analysis event we write below. Every value is
# stringified on the way out; metadata is a mapping<string, string> in the
# GenomeAnnotation spec, so a consumer never has to care whether a count arrived
# as a number or a string, and JSON round-trips it identically either way.
#
my %meta = ( input_contigs => scalar @{$genome_in->{contigs}} );

my $reference_name = $opt->reference;
my $reference_source;
my $taxon;

if ($reference_name)
{
    $reference_source = 'option';
}
else
{
    $taxon = $opt->taxon // $genome_in->taxonomy_id;
    if ($taxon)
    {
	my $api = P3DataAPI->new;
	$reference_name = Bio::BVBRC::ViralAnnotation::VigorTaxonMap::find_vigor_reference($taxon, $api);
	if ($reference_name)
	{
	    $reference_source = 'taxon_map';
	}
	else
	{
	    warn "No reference found for taxon $taxon\n";
	}
    }
}

$meta{taxon}            = $taxon            if $taxon;
$meta{reference_name}   = $reference_name   if $reference_name;
$meta{reference_source} = $reference_source if $reference_source;

if (!$reference_name)
{
    warn "No reference found\n";

    #
    # Record the run even though we are not going to do anything. Without an event
    # here, a no-reference run is indistinguishable from the stage never having been
    # scheduled at all -- which is exactly the distinction a pipeline condition
    # gating a downstream stage needs to make.
    #
    $meta{status} = 'no_reference';
    my($event) = make_event($genome_in, []);
    finish_event($event, \%meta, 0);

    $genome_in->destroy_to_file($opt->output);
    exit 0;
}

#
# Invoke vigor4 to annotate viral genome.
#
# Write contigs as fasta file.
# Invoke vigor4, using the reference that was passed in as the reference_name parameter;
# We parse the .pep file that is generated. It contains two feature types.
# This is a CDS:
# >NC_045512.1 location=266..13468,13471..21555 codon_start=1 gene="orf1ab" ref_db="covid19" ref_id="YP_009724389.1"
# This is a mature peptide:
# >NC_045512.1.1 mat_peptide location=266..805 gene="orf1ab" product="leader protein" ref_db="covid19_orf1ab_mp" ref_id="YP_009725297.1"
#
# We create features of type CDS and mat_peptide.
#

my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();

my $ref_dir = vigor_reference_db_directory;
$ref_dir ne '' or die "Vigor reference directory not configured";
-d $ref_dir or die "Vigor reference directory '$ref_dir' not found";

my @vigor_params = ("-i", $sequences_file,
		    "--reference-database-path", $ref_dir,
		    "-d", $reference_name,
		    "-o", "$here/vigor_out");

print STDERR Dumper(\@vigor_params);
my $ok = run(["vigor4", @vigor_params],
	     init => sub {
		 chdir $tempdir;
		 $ENV{JAVA_OPTS} = "-XX:ActiveProcessorCount=" . $opt->threads;
	     },
	     ">", "$here/vigor4.stdout.txt",
	     "2>", "$here/vigor4.stderr.txt");
my $vigor_rc = $?;

#
# Report the exit code, not the raw wait status -- $? is 768 for an exit(3), which
# is needlessly confusing to whoever reads this metadata.
#
$meta{exit_code} = $vigor_rc >> 8;
$meta{signal}    = $vigor_rc & 127 if $vigor_rc & 127;

if (!$ok)
{
    print STDERR "Vigor run failed with rc=$vigor_rc. Stdout:\n";
    copy("$here/vigor4.stdout.txt", \*STDERR);
    print STDERR "Stderr:\n";
    copy("$here/vigor4.stderr.txt", \*STDERR);
}
    
my($event, $event_id) = make_event($genome_in, \@vigor_params);

#
# "Success" is the question a caller actually has: did this run produce an
# annotation? Not "did the process exit 0" -- vigor can exit non-zero and still
# leave a parseable .pep -- and not "does the GTO have features", which is true of
# any GTO that arrived carrying GenBank features.
#
$meta{status} = $ok ? 'ok' : 'vigor_failed';

#
# Always report all three counts, including the zeroes. A reference database for a
# non-polyprotein virus (influenza, rotavirus, RSV, the Bunyavirales genus dbs)
# legitimately yields no mature peptides, and a condition downstream has to be able
# to tell that from "vigor4 never reported". An absent key cannot say it.
#
my %counts = (CDS => 0, mat_peptide => 0, pseudogene => 0);

#
# Parse the generated peptide file. We collect the CDS and mature_peptides, then
# add features so that we can register the counts.
#
if (open(my $pep_fh, "<", "$here/vigor_out.pep"))
{
    my %features;
    while (my($id, $def, $seq) = read_next_fasta_seq($pep_fh))
    {
	my $fq = { truncated_begin => 0, truncated_end => 0 };
	
	my $type;
	my $ctg;
	if ($def =~ s/^mat_peptide\s+//)
	{
	    ($ctg) = $id =~ /^(.*)\.[^.]+\.[^.]+$/;
	    $type = 'mat_peptide';
	}
	elsif ($def =~ s/^pseudogene\s+//)
	{
	    ($ctg) = $id =~ /^(.*)\.[^.]+$/;
	    $type = 'pseudogene';
	}
	else
	{
	    ($ctg) = $id =~ /^(.*)\.[^.]+$/;
	    $type = 'CDS';
	}
	
	if (!$ctg)
	{
	    print STDERR "Falling back to prefix of id for contig name from $id\n";
	    ($ctg) = $id =~ /^(.*?)\./;
	}
	
	my $feature = {
	    quality => $fq,
	    type => $type,
	    contig => $ctg,
	    aa_sequence => $seq,
	};
	push(@{$features{$type}}, $feature);
	
	while ($def =~ /([^=]+)=((\"([^\"]+)\")|([^\"\s]+))\s*/mg) 
	{
	    my $key = $1;
	    my $val = $4 ? $4 : $5;
	    
	    my @loc;
	    # print "key=$key val=$val\n";

	    if ($key eq 'location')
	    {
		$feature->{genbank_feature} = { genbank_type => $type, genbank_location  => $val, values => {}};
		
		# location=266..13468,13471..21555
		for my $ent (split(/,/, $val))
		{
		    if (my($s_frag, $s, $e_frag, $e) = $ent =~ /^(<?)(\d+)\.\.(>?)(\d+)$/)
		    {
			$fq->{truncated_begin} = 1 if $s_frag;
			$fq->{truncated_end} = 1 if $e_frag;
			
			my $len = abs($s - $e) + 1;
			my $strand = $s < $e ? '+' : '-';
			push(@loc, [$ctg, $s, $strand, $len]);
		    }
		    else
		    {
			die "error parsing location '$ent'\n";
		    }
		}
		$feature->{location} = \@loc;
	    }
	    else
	    {
		$feature->{$key} = $val;
	    }
	}
	$feature->{product} //= $feature->{gene};
	if (!$feature->{location})
	{
	    warn "No location for feature $def " . Dumper($feature);
	}
    }
    # print Dumper(BEFORE => \%features);

    #
    # If we have features, remove any existing CDS or mat_peptide features.
    #

    if (%features && $opt->remove_existing)
    {
	my @to_del = $genome_in->fids_of_type('CDS', 'mat_peptide', 'pseudogene');
	print STDERR "Delete @to_del\n";
	$meta{removed_existing} = scalar @to_del;
	$genome_in->delete_feature($_) foreach @to_del;
    }
    # print Dumper(AFTER => $genome_in);
    
    for my $type (keys %features)
    {
	my $feats = $features{$type};
	my $n = @$feats;
	my $id_type = $type;

	$counts{$type} = $n;
	
	for my $feature (@$feats)
	{
	    my $p = {
		-id	     => $genome_in->new_feature_id($id_type),
		-type 	     => $type,
		-location 	     => $feature->{location},
		-analysis_event_id 	     => $event_id,
		-annotator => 'vigor4',
		-protein_translation => $feature->{aa_sequence},
		-alias_pairs => [[gene => $feature->{gene}]],
		-function => $feature->{product},
		-quality_measure => $feature->{quality},
		-genbank_feature => $feature->{genbank_feature},
	    };
	    #die Dumper($p);
	    
	    $genome_in->add_feature($p);
	}
    }
    
}
else
{
    warn "Could not read $here/vigor_out.pep\n";
    $meta{status} = 'no_pep_file';
}

my $features_called = 0;
for my $type (sort keys %counts)
{
    $meta{"${type}_called"} = $counts{$type};
    $features_called += $counts{$type};
}
$meta{features_called} = $features_called;

finish_event($event, \%meta, $features_called);

$genome_in->destroy_to_file($opt->output);

#
# Create the analysis event and register it with the GTO.
#
# add_analysis_event stores the hashref we hand it rather than a copy, so the event
# is registered up front -- the features we add need its id -- and success and
# metadata are filled in by finish_event() below, once there is something to
# report. The change is reflected in the GTO that gets written.
#
sub make_event
{
    my($gto, $params) = @_;

    my $event = {
	tool_name => "vigor4",
	execution_time => scalar gettimeofday,
	parameters => $params,
	hostname => $hostname,
    };

    my $event_id = $gto->add_analysis_event($event);

    return($event, $event_id);
}

sub finish_event
{
    my($event, $meta, $success) = @_;

    $event->{success} = $success ? 1 : 0;
    $event->{metadata} = { map { $_ => "$meta->{$_}" } grep { defined $meta->{$_} } keys %$meta };
}

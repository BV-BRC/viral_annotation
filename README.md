# Viral Annotation Tools

## Overview

This module contains BV-BRC genome annotation service wrapper scripts for two viral annotation tools: VIGOR4 and the VIPR Mature Peptide annotation tool.

## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md). The scripts in this module are invoked by the annnotation pipeline found [in this module](https://github.com/TheSEED/genome_annotation).

## VIGOR4

VIGOR4 (Viral Genome ORF Reader) is a Java application to predict protein sequences encoded in viral genomes.
VIGOR4 determines the protein coding sequences by sequence similarity searching against curated viral protein databases.

It utilizes the [VIGOR_DB](https://github.com/JCVenterInstitute/VIGOR_DB) project which maintains reference databases for the viruses enumerated in the [vigor-taxon-map](vigor-taxon-map.txt) file.

The first reference below describes the overall implementation of the VIGOR tool.

### References

Wang S, Sundaram JP, Spiro D. VIGOR, an annotation program for small viral genomes. BMC Bioinformatics. 2010 Sep 7;11:451. doi: 10.1186/1471-2105-11-451. PMID: 20822531; PMCID: PMC2942859.

Wang S, Sundaram JP, Stockwell TB. VIGOR extended to annotate genomes for additional 12 different viruses. Nucleic Acids Res. 2012 Jul;40(Web Server issue):W186-92. doi: 10.1093/nar/gks528. Epub 2012 Jun 4. PMID: 22669909; PMCID: PMC3394299.

## Mature Peptide Annotation

The following discussion is extracted and lightly edited from [the module providing the underlying analysis code](https://github.com/BV-BRC/vipr_mat_peptide/blob/master/vipr_mat_peptide_readme.txt).

Some viruses code all propteins in one long coding sequence (CDS), which produces a polyprotein after translation.
The polyprotein is then processed to generate final products, the mature peptides, by proteinase. These mat_peptides are the real functional proteins.
Many viral genomes in genbank lack the mat_peptide annotation, which is essential to study the sequence, structure, and function of a protein, as mat_peptide is the functional unit.
Other genomes might have mat_peptide annotation that doesn't confirm with latest development

### Using reference sequences to annotate mat_peptides

Many viral species have reference sequences in Genbank with mat_peptide annotation, which can be used as standard to annotate other genomes of the same species.
These reference sequences have been created by expert biologist as well as NCBI staff.
Some have undergone review, while others are provisional.
Most have complete sequence in full length

### Annotation Process

With a normal genome, one would do following:

1. find an appropriate reference sequence,
2. find the corresponding CDS,
3. align the CDS, then
4. find the start/end of each mat_peptide in the target,
5. define the mat_peptides in the target, and finally
6. save the annotated mat_peptide.

### Workflow of the Script

1. find an appropriate reference sequence:
   - First determine the taxon id of the target genome, then look for a defined reference sequence for that strain.
   - If no reference sequence for that strain, go up to the species level and look for any defined reference sequence for the species
2. find the corresponding CDS:
   - There might be >1 CDSs in either target or reference sequence. To find the corresponding CDSs, run blast with a target CDS and a reference sequence CDS. Anything with <50% conserved residues is ignored
3. align the CDS:
   - Run global alignment with CLUSTALW (MUSCLE was used before Oct 2012).
   - Blast produces local alignment, which might leave out any overhang, or mis-matched terminals
     -Also considered hidden Markov model, but this is more of a scoring method than a alignment method
4. find the start/end of each mat_peptide in the target:
   - Based on the alignment of target CDS and reference sequence CDS, find the start/end of each mat_peptide in the target genome, taking into consideration following factors: different start of the polyprotein, terminal gaps in MSA, and internal gaps in MSA, any skip/buckle in either reference sequence/target CDS. Then convert the AA range back to DNA range
5. define the mat_peptides in the target:
   - For each mat_peptide, identify with genome accession (with version), CDS id, reference sequence accession, DNA location in genome, AA location in polyprotein, gene symbol, and product name
6. save the annotated mat_peptide:
   - Save the annotation in fasta for database loading. May also save updated genbank to file for other purpose

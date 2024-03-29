TOP_DIR = ../..
include $(TOP_DIR)/tools/Makefile.common

DEPLOY_RUNTIME ?= /kb/runtime
TARGET ?= /kb/deployment

APP_SERVICE = app_service

SRC_PERL = $(wildcard scripts/*.pl)
BIN_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_PERL))))
DEPLOY_PERL = $(addprefix $(TARGET)/bin/,$(basename $(notdir $(SRC_PERL))))

SRC_SERVICE_PERL = $(wildcard service-scripts/*.pl)
BIN_SERVICE_PERL = $(addprefix $(BIN_DIR)/,$(basename $(notdir $(SRC_SERVICE_PERL))))
DEPLOY_SERVICE_PERL = $(addprefix $(SERVICE_DIR)/bin/,$(basename $(notdir $(SRC_SERVICE_PERL))))

CLIENT_TESTS = $(wildcard t/client-tests/*.t)
SERVER_TESTS = $(wildcard t/server-tests/*.t)
PROD_TESTS = $(wildcard t/prod-tests/*.t)

STARMAN_WORKERS = 8
STARMAN_MAX_REQUESTS = 100

VIGOR_REFERENCE_DB_DIRECTORY = /disks/patric-common/runtime/vigor-4.1.20190809-121720-7fa683e/VIGOR_DB

TPAGE_ARGS = --define kb_top=$(TARGET) --define kb_runtime=$(DEPLOY_RUNTIME) --define kb_service_name=$(SERVICE) \
	--define kb_service_port=$(SERVICE_PORT) --define kb_service_dir=$(SERVICE_DIR) \
	--define kb_sphinx_port=$(SPHINX_PORT) --define kb_sphinx_host=$(SPHINX_HOST) \
	--define kb_starman_workers=$(STARMAN_WORKERS) \
	--define kb_starman_max_requests=$(STARMAN_MAX_REQUESTS) \
	--define vigor_reference_db_directory=$(VIGOR_REFERENCE_DB_DIRECTORY) \
	--define viral_family_db=$(VIRAL_FAMILY_DB)

all: build-libs bin 

build-libs:
	mkdir -p lib/Bio/BVBRC/ViralAnnotation
	$(TPAGE) $(TPAGE_BUILD_ARGS) $(TPAGE_ARGS) Config.pm.tt > lib/Bio/BVBRC/ViralAnnotation/Config.pm
	$(PERL) build-taxon-map.pl  vigor-taxon-map.txt VigorTaxonMap.pm.tt lib/Bio/BVBRC/ViralAnnotation/VigorTaxonMap.pm

bin: $(BIN_PERL) $(BIN_SERVICE_PERL)

deploy: deploy-all
deploy-all: deploy-client 
deploy-client: build-libs deploy-libs deploy-scripts deploy-docs

deploy-service: build-libs deploy-libs deploy-scripts deploy-service-scripts deploy-specs

deploy-dir:
	if [ ! -d $(SERVICE_DIR) ] ; then mkdir $(SERVICE_DIR) ; fi
	if [ ! -d $(SERVICE_DIR)/bin ] ; then mkdir $(SERVICE_DIR)/bin ; fi

deploy-docs: 


clean:



include $(TOP_DIR)/tools/Makefile.common.rules

### Compiler ###


GCC = $(shell which gcc)
GPP = $(shell which g++)
ICC = $(shell which icc)
CXX = $(GPP)


### CS-SUITE version ###


VERSION = 2.2.3
OUTFILE = $@
ifdef SHOW_VERSION
	OUTFILE = $@-$(VERSION)
endif


### Flags ###


FLAGS =						# Global flags
RELEASE_FLAGS =		# Flags for debug targets
DEBUG_FLAGS =			# Flags for release targets
GPROF_FLAGS = -g -pg -fno-omit-frame-pointer -fno-inline-functions -DNDEBUG
GTEST_FLAGS = -I$(GTEST_DIR) -I$(GTEST_DIR)/include -I/usr/include -DGTEST_USE_OWN_TR1_TUPLE=1 -DGTEST_HAS_PTHREAD=0

FLAGS = -Wall -Wno-deprecated
ifeq ($(BITS), 32)
	FLAGS += -m32
else
	FLAGS += -m64
endif

RELEASE_FLAGS = -O3 -DNDEBUG -fopenmp -DOPENMP $(FLAGS)
DEBUG_FLAGS = -g3 $(FLAGS)

ifeq ($(OPENMP),1)
	DEBUG_FLAGS += -fopenmp -DOPENMP
endif


### Definitions ###


DEFS =	-fno-strict-aliasing \
		  -DVERSION="\"$(VERSION)\"" \
     	-DBUILD_HOST="\"`hostname`\"" \
     	-DBUILD_TIME="\"`date`\"" \
     	-DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\""
ifdef LOGGING
	DEFS += -DLOGGING -DLOG_MAX_LEVEL=$(LOGGING) 
endif

RELEASE_DEFS = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS)\"" $(DEFS)
DEBUG_DEFS = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS)\"" $(DEFS)


### Files ###


BUILD_ROOT = ../build
ifeq ($(MODE),test)
	BUILD_DIR = $(BUILD_ROOT)/test
else ifeq  ($(MODE),debug)
	BUILD_DIR = $(BUILD_ROOT)/debug
else ifeq  ($(MODE),prof)
	BUILD_DIR = $(BUILD_ROOT)/prof
else
	BUILD_DIR = $(BUILD_ROOT)
endif
$(shell mkdir -p $(BUILD_DIR))
BIN_DIR = ../bin
$(shell mkdir -p $(BIN_DIR))
LIB_DIR = ../lib
GTEST_DIR = ../lib/gtest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRC = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
EXT = cc
HEADERS = $(wildcard *.h)
DEPS =
BASE_DEPS = application aa dna log getopt_pp blosum_matrix assert_helpers utils
GTEST_DEPS = gtest-all
ifeq ($(MODE),test)
	OBJECTS = $(patsubst %, $(BUILD_DIR)/%.o, $(BASE_DEPS) $(DEPS) $(GTEST_DEPS))
else
	OBJECTS = $(patsubst %, $(BUILD_DIR)/%.o, $(BASE_DEPS) $(DEPS))
endif
INC = -I$$HOME/include
LIBS = -lm


### Targets ###


TARGETS = csblast cstrainset cssgd csbuild csviz cstranslate cscons \
					cscp_neff cstrainset_neff csclust
BINS = $(TARGETS:%=$(BIN_DIR)/%)

ifeq ($(MODE),test)
	PARAMS = $(DEBUG_FLAGS) $(DEBUG_DEFS) $(INC) $(GTEST_FLAGS)
else ifeq ($(MODE),debug)
	PARAMS = $(DEBUG_FLAGS) $(DEBUG_DEFS) $(INC)
else ifeq ($(MODE),prof)
	PARAMS = $(RELEASE_FLAGS) $(RELEASE_DEFS) $(INC) $(GPROF_FLAGS)
else
	PARAMS = $(RELEASE_FLAGS) $(RELEASE_DEFS) $(INC)
endif


### Global targets ###


all: $(TARGETS)

clean_build:
	find $(BUILD_DIR) -name '*.o' -exec rm -f '{}' +

clean_bin:
	rm -f $(BINS)

clean: clean_build clean_bin

.PHONY: all all_test clean clean_build clean_bin

$(BUILD_DIR)/%.o: %.$(EXT) $(HEADERS)
	$(CXX) $(PARAMS) -c -o $@ $<


### GTEST targets


$(BUILD_DIR)/gtest-all.o : $(GTEST_SRC)
	$(CXX) $(PARAMS) -c -o $@ $(GTEST_DIR)/src/gtest-all.cc

$(BUILD_DIR)/gtest_main.o : $(GTEST_SRC)
	$(CXX) $(PARAMS) -c -o $@ $(GTEST_DIR)/src/gtest_main.cc


### csblast ###


DEPS = csblast_app csblast_iteration csblast blast_hits
csblast: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cstrainset ###


DEPS = cstrainset_app
cstrainset: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cssgd ###


DEPS = cssgd_app
cssgd: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cshmc ###


DEPS = cshmc_app
cshmc: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csbuild ###


DEPS = csbuild_app
csbuild: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cscons ###


DEPS = cscons_app
cscons: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csviz ###


DEPS = csviz_app pdf_writer as
csviz: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csclust ###


DEPS = csclust_app pdf_writer as
csclust: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cstranslate ###


DEPS = cstranslate_app as
cstranslate: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cscp_neff ###


DEPS = cscp_neff_app as
cscp_neff: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cstrainset_neff ###


DEPS = cstrainset_neff_app as
cstrainset_neff: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csconvert ###


DEPS = csconvert_app as
csconvert: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### Test targets ###


DEPS = test
test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = sgd_test
sgd_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = crf_test
crf_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = count_profile_test
count_profile_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = profile_test
profile_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = library_pseudocounts_test
library_pseudocounts_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = pseudocounts_test
pseudocounts_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = alignment_test blast_hits
alignment_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

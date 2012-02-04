# Author	: Angermueller Christof
# Date		: 10-09-09


### Compiler ###


GCC = $(shell which gcc)
GPP = $(shell which g++)
ICC = $(shell which icc)
CXX = $(GPP)


### CS-SUITE version ###


VERSION = 2.2.0
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
ifdef OPENMP
	FLAGS += -fopenmp -DOPENMP
endif

RELEASE_FLAGS = -O3 -DNDEBUG -fopenmp -DOPENMP $(FLAGS)
DEBUG_FLAGS = -g3 $(FLAGS)


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


BUILD_DIR = ../build
BUILD_DIR_TEST = $(BUILD_DIR)/test
BUILD_DIR_PROF = $(BUILD_DIR)/prof
BUILD_DIR_DEBUG = $(BUILD_DIR)/debug
$(shell mkdir -p $(BUILD_DIR))
$(shell mkdir -p $(BUILD_DIR_TEST))
$(shell mkdir -p $(BUILD_DIR_DEBUG))
$(shell mkdir -p $(BUILD_DIR_PROF))
BIN_DIR = ../bin
LIB_DIR = ../lib
GTEST_DIR = ../lib/gtest
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_SRC = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)
EXT = cc
HEADERS = $(wildcard *.h)
DEPS =
BASE_DEPS = application aa dna log getopt_pp blosum_matrix assert_helpers utils
GTEST_DEPS = gtest_main gtest-all
OBJECTS = $(patsubst %, $(BUILD_DIR)/%.o, $(BASE_DEPS) $(DEPS))
OBJECTS_TEST = $(patsubst %, $(BUILD_DIR_TEST)/%.o, $(BASE_DEPS) $(DEPS) $(GTEST_DEPS))
OBJECTS_PROF = $(patsubst %, $(BUILD_DIR_PROF)/%.o, $(BASE_DEPS) $(DEPS))
OBJECTS_DEBUG = $(patsubst %, $(BUILD_DIR_DEBUG)/%.o, $(BASE_DEPS) $(DEPS))
INC = -I$$HOME/lib/sparsehash/include
LIBS = -lm


### Targets ###


TARGETS = cslrsgd csblast cstrainset cssgd csbuild csviz csalign csclust cstranslate \
					cscp_neff cstrainset_neff cstrainprofiles
TEST_TARGETS = lr_params lr_func lr_sgd lr_pseudocounts training_profile
TARGETS := $(TARGETS) $(TARGETS:%=%_debug) $(TARGETS:%=%_prof)
TEST_TARGETS := $(TEST_TARGETS:%=%_test)
BINS = $(TARGETS:%=$(BIN_DIR)/%) $(TEST_TARGETS:%=$(BIN_DIR)/%)

PARAMS = $(RELEASE_FLAGS) $(RELEASE_DEFS) $(INC)
%_test: PARAMS = $(DEBUG_FLAGS) $(DEBUG_DEFS) $(INC) $(GTEST_FLAGS)
%_debug: PARAMS = $(DEBUG_FLAGS) $(DEBUG_DEFS) $(INC)
%_prof: PARAMS = $(RELEASE_FLAGS) $(RELEASE_DEFS) $(INC) $(GPROF_FLAGS)


### Global targets ###


all: $(TARGETS)

all_test: $(TEST_TARGETS)

clean_build:
	find $(BUILD_DIR) -name '*.o' -exec rm -f '{}' +

clean_build_test:
	rm $(BUILD_DIR_TEST)/*.o

clean_build_prof:
	rm $(BUILD_DIR_PROF)/*.o

clean_build_debug:
	rm $(BUILD_DIR_DEBUG)/*.o

clean_bin:
	rm -f $(BINS)

clean: clean_build clean_build_test clean_build_prof clean_build_debug cclean_bin

.PHONY: all all_test clean clean_build clean_bin

$(BUILD_DIR)/%.o: %.$(EXT) $(HEADERS)
	$(CXX) $(PARAMS) -c -o $@ $<

$(BUILD_DIR_TEST)/%.o: %.$(EXT) $(HEADERS)
	$(CXX) $(PARAMS) -c -o $@ $<

$(BUILD_DIR_PROF)/%.o: %.$(EXT) $(HEADERS)
	$(CXX) $(PARAMS) -c -o $@ $<

$(BUILD_DIR_DEBUG)/%.o: %.$(EXT) $(HEADERS)
	$(CXX) $(PARAMS) -c -o $@ $<


### GTEST targets


$(BUILD_DIR_TEST)/gtest-all.o : $(GTEST_SRC)
	$(CXX) $(PARAMS) -c -o $@ $(GTEST_DIR)/src/gtest-all.cc

$(BUILD_DIR_TEST)/gtest_main.o : $(GTEST_SRC)
	$(CXX) $(PARAMS) -c -o $@ $(GTEST_DIR)/src/gtest_main.cc


### csblast ###


DEPS = csblast_app csblast_iteration csblast blast_hits
csblast: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csblast_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csblast_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cstrainset ###


DEPS = cstrainset_app
cstrainset: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cstrainset_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cstrainset_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cssgd ###


DEPS = cssgd_app
cssgd: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cssgd_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cssgd_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cslrsgd ###


DEPS = cslrsgd_app
cslrsgd: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cslrsgd_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cslrsgd_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cshmc ###


DEPS = cshmc_app
cshmc: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cshmc_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cshmc_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csbuild ###


DEPS = csbuild_app
csbuild: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csbuild_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csbuild_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csviz ###


DEPS = csviz_app pdf_writer as
csviz: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csviz_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csviz_prof: $(OBJECTS_RPOF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csalign ###


DEPS = csalign pdf_writer sparse_profile
csalign: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csalign_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csalign_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### csclust ###


DEPS = csclust_app pdf_writer as
csclust: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csclust_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

csclust_prof: $(OBJECTS_PROF)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)


### cstranslate ###


DEPS = cstranslate_app as
cstranslate: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cstranslate_debug: $(OBJECTS_DEBUG)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

cstranslate_prof: $(OBJECTS_PROF)
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

DEPS = lr_params_test
lr_params_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = lr_func_test
lr_func_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = lr_sgd_test
lr_sgd_test: $(OBJECTS_TEST)
	$(CXX) $(PARAMS) -o $(BIN_DIR)/$(OUTFILE) $+ $(LIBS)

DEPS = lr_pseudocounts_test
lr_pseudocounts_test: $(OBJECTS_TEST)
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

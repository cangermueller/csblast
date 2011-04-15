# Author	: Angermueller Christof
# Date		: 10-09-09


### Compiler ###


GCC_PREFIX = $(shell dirname `which gcc`)
GCC_SUFFIX =
CC = $(GCC_PREFIX)/gcc$(GCC_SUFFIX)
CPP = $(GCC_PREFIX)/g++$(GCC_SUFFIX)
CXX = $(CPP)


### Flags ###


FLAGS =				# Global flags
RELEASE_FLAGS =		# Flags for debug targets
DEBUG_FLAGS =		# Flags for release targets
GPROF_FLAGS = -g -pg -fno-omit-frame-pointer -fno-inline-functions -DNDEBUG
GTEST_FLAGS = $(shell gtest-config --libs --cppflags)
LOG_MAX_LEVEL = 0

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
DEBUG_FLAGS = -O1 -g3 $(FLAGS)


### Definitions ###


DEFS =	-fno-strict-aliasing \
     	-DBUILD_HOST="\"`hostname`\"" \
     	-DBUILD_TIME="\"`date`\"" \
     	-DCOMPILER_VERSION="\"`$(CXX) -v 2>&1 | tail -1`\""
ifdef LOGGING
	DEFS += -DLOGGING -DLOG_MAX_LEVEL=$(LOG_MAX_LEVEL) 
endif

RELEASE_DEFS = -DCOMPILER_OPTIONS="\"$(RELEASE_FLAGS)\"" $(DEFS)
DEBUG_DEFS = -DCOMPILER_OPTIONS="\"$(DEBUG_FLAGS)\"" $(DEFS)


### Files ###


BUILDDIR = ../build
ifdef BUILD
	BUILDDIR := $(BUILDDIR)/$(BUILD)
    $(shell mkdir -p $(BUILDDIR))
endif
BINDIR = ../bin
EXT = cc
HEADERS = $(wildcard *.h)
DEPS =
BASE_DEPS = application aa dna log getopt_pp blosum_matrix assert_helpers utils
OBJECTS = $(patsubst %, $(BUILDDIR)/%.o, $(BASE_DEPS) $(DEPS))
INC =
LIBS = -lm

### Targets ###


TARGETS = cslrsgd csblast cstrainset cssgd csbuild csviz csalign csclust cstranslate \
					count_profile_neff trainset_neff
TEST_TARGETS = lr_params lr_func lr_sgd lr_pseudocounts training_profile
TARGETS := $(TARGETS) $(TARGETS:%=%_debug) $(TARGETS:%=%_prof)
TEST_TARGETS := $(TEST_TARGETS:%=%_test)
BINS = $(TARGETS:%=$(BINDIR)/%) $(TEST_TARGETS:%=$(BINDIR)/%)

PARAMS = $(RELEASE_FLAGS) $(RELEASE_DEFS) $(INC)
%_test: PARAMS = $(DEBUG_FLAGS) $(DEBUG_DEFS) $(INC) $(GTEST_FLAGS)
%_debug: PARAMS = $(DEBUG_FLAGS) $(DEBUG_DEFS) $(INC)
%_prof: PARAMS = $(RELEASE_FLAGS) $(RELEASE_DEFS) $(INC) $(GPROF_FLAGS)


### Global targets ###


all: $(TARGETS)

all_test: $(TEST_TARGETS)

clean_build:
	find $(BUILDDIR) -name '*.o' -exec rm -f '{}' +

clean_bin:
	rm -f $(BINS)

clean: clean_build clean_bin

.PHONY: all all_test clean clean_build clean_bin

$(BUILDDIR)/%.o: %.$(EXT) $(HEADERS)
	$(CXX) $(PARAMS) -c -o $@ $<


### Test targets ###


DEPS = cs_test 
cs_test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

DEPS = lr_params_test
lr_params_test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

DEPS = lr_func_test
lr_func_test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

DEPS = lr_sgd_test
lr_sgd_test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

DEPS = lr_pseudocounts_test
lr_pseudocounts_test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

DEPS = training_profile_test
training_profile_test: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### cslrsgd ###


DEPS = cslrsgd_app
cslrsgd: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cslrsgd_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cslrsgd_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### csblast ###


DEPS = csblast_app csblast_iteration csblast blast_hits
csblast: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csblast_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csblast_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### cstrainset ###


DEPS = cstrainset_app
cstrainset: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cstrainset_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cstrainset_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### cssgd ###


DEPS = cssgd_app
cssgd: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cssgd_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cssgd_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### csbuild ###


DEPS = csbuild_app
csbuild: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csbuild_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csbuild_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### csviz ###


DEPS = csviz_app pdf_writer as
csviz: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csviz_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csviz_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### csalign ###


DEPS = csalign pdf_writer sparse_profile
csalign: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csalign_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csalign_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### csclust ###


DEPS = csclust_app pdf_writer as
csclust: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csclust_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

csclust_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### cstranslate ###


DEPS = cstranslate_app as
cstranslate: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cstranslate_debug: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

cstranslate_prof: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### count_profile_neff ###


DEPS = count_profile_neff_app as
count_profile_neff: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)


### trainset_neff ###


DEPS = trainset_neff_app as
trainset_neff: $(OBJECTS)
	$(CXX) $(PARAMS) -o $(BINDIR)/$@ $+ $(LIBS)

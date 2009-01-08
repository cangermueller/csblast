TARGETS  :=
CXXFLAGS := -O3 -g -W -Wall -Wno-long-long -Wno-char-subscripts -fno-strict-aliasing -I.
CXX      := g++
LIBS     :=
EXT      := cpp
BUILDDIR := build

override BUILDDIR := $(strip $(BUILDDIR))
SOURCES  := $(wildcard *.$(EXT))
OBJECTS  := $(patsubst %.$(EXT), $(BUILDDIR)/%.o, $(SOURCES))
DEPS     := $(patsubst %.$(EXT), $(BUILDDIR)/%.dep, $(SOURCES))
TESTS    := $(wildcard *test.h)

# If you want to use python, specify USE_PYTHON=1 on the command line
ifdef USE_PYTHON
	TESTGEN = ../cxxtestgen.py
else
	TESTGEN = ../cxxtestgen.pl
endif

.PHONY: all
all: $(TARGETS)

# xxx: $(OBJECTS) $(DEPS)
# 	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(LIBS)

ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

$(OBJECTS): $(BUILDDIR)/%.o: %.$(EXT) $(BUILDDIR)/%.dep $(BUILDDIR)/.tag
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(DEPS): $(BUILDDIR)/%.dep: %.$(EXT) $(BUILDDIR)/.tag
	$(CXX) $(CXXFLAGS) -MM $< -MT $@ -MT $(<:.$(EXT)=.o) -o $@

%.tag:
	mkdir -p $(dir $(@))
	touch $@

run_tests: runner
	./runner

runner: runner.cpp $(OBJECTS) $(DEPS)
	$(CXX) $(CXXFLAGS) -o $@ runner.cpp $(OBJECTS) $(LIBS)

runner.cpp: $(TESTS)
	$(TESTGEN) -o $@ --error-printer $(TESTS)

.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(TARGETS) runner runner.cpp

CXX = g++
CXXFLAGS = -O3 -Wall -Wno-char-subscripts -fno-strict-aliasing
CPP = $(CXX) $(CXXFLAGS)

SOURCES = test.cpp Sequence_profile.cpp Amino_acid_alphabet.cpp Nucleic_acid_alphabet.cpp Sequence_alphabet.cpp 
HEADERS = Sequence_profile.h Amino_acid_alphabet.h Nucleic_acid_alphabet.h Sequence_alphabet.h

TARGETS = test

all: $(TARGETS)

test: $(SOURCES) $(HEADERS)
	$(CPP) -o ../bin/test $(SOURCES)

clean: 	
	-rm -f *.o ../bin/*

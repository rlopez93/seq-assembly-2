# CXX=g++
CXX_FLAGS+=-Wall -Wextra -Wpedantic -std=c++1y -O3

all: assemble

clean:
	rm -rf assemble khmer/*.o

assemble: khmer/kmer_hash.o khmer/MurmurHash3.o assemble.cpp
	$(CXX) $(CXX_FLAGS) assemble.cpp khmer/kmer_hash.o khmer/MurmurHash3.o -o assemble

khmer/kmer_hash.o: khmer/kmer_hash.cc khmer/kmer_hash.hh
	$(CXX) $(CXX_FLAGS) khmer/kmer_hash.cc -c -o khmer/kmer_hash.o

khmer/MurmurHash3.o: khmer/MurmurHash3.cc khmer/MurmurHash3.h
	$(CXX) $(CXX_FLAGS) khmer/MurmurHash3.cc -c -o khmer/MurmurHash3.o

#! /usr/bin/env python2
"""
--- assemble.py ---

Usage: python assemble.py <FAST[AQ]> <ksize>

Receives DNA read data as input in FASTA or FASTQ format
and a k-mer size, and prints the assembled sequences
to stdout in FASTA format.

Example:

$ python assemble.py sample.fasta 18 > assembled.fasta
"""
from __future__ import print_function
from euler import find_eulerian_path
import sys
import screed
import networkx as nx

if len(sys.argv) != 3:
    print(__doc__)
    sys.exit(1)

# input file containing reads
# input may be in FASTA or FASTQ format
infilename = sys.argv[1]

# kmer size to use in building De Bruijn Graph
k = int(sys.argv[2])

# De Bruijn Graph of k-mers
# using a NetworkX DiGraph
DG = nx.DiGraph()

# iterate over reads in input file using screed
for record in screed.open(infilename):
    seq = record.sequence  # get current read

    # iterate over all k-mers in seq,
    # and add them to the graph
    for i in xrange(len(seq)-k-1):
        kmer_a = seq[i:i+k]
        kmer_b = seq[i+1:i+k+1]
        DG.add_edge(kmer_a, kmer_b)


# function for printing sequence
# with 60 bases per row
def print_seq(a, seq):
    i = 0
    for base in a:
        print(base, end='')
        if (i+1) % 60 == 0:
            print()
        i += 1

    for base in seq:
        print(base, end='')
        if (i+1) % 60 == 0:
            print()
        i += 1
    print()

# generator for weakly connected component subgraphs of DG
# each represents a possible contig in the original sequence
wc_subgraphs = nx.weakly_connected_component_subgraphs(DG)

for i, subgraph in enumerate(wc_subgraphs):

    # try block will throw if subgraph
    # does not contain an eulerian path
    try:
        path = find_eulerian_path(subgraph)

        # find first two vertices in path
        a, b = path.next()
        first = a + b[-1]

        # get rest of path using a generator expression
        subseq = (c[-1] for _, c in path)

        # print sequence in fasta format
        print(">{}".format(i))
        print_seq(first, subseq)

    except nx.NetworkXError:
        # subgraph does not contain an eulerian path
        # so do nothing
        pass

#! /usr/bin/env python
"""
Usage: python assemble.py <FAST[AQ]> <ksize>
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
        path.next()

    except nx.NetworkXError:
        # subgraph does not contain an eulerian path

        # what follows is code that tries to extract
        # possible sequences from the non-eulerian subgraph

        # to do that, I set all nodes with in degree 0
        # as possible sources and all nodes with out degree 0
        # as possible targets.

        # I then find a shortest path between all sources
        # and all targets, if it exists, and print those
        # paths as the contigs

        # WARNING - HACK APPROACHING
        # uncomment if you dare

        sources = []
        targets = []

        for node in subgraph:
            ind = subgraph.in_degree(node)
            outd = subgraph.out_degree(node)
            if ind == 0:
                sources.append(node)
            if outd == 0:
                targets.append(node)

        for j, source in enumerate(sources):
            for k, target in enumerate(targets):
                # try will throw if there is no path between source and target
                try:
                    sp = nx.shortest_path(subgraph, source, target)

                    first = sp[0]
                    subseq = (v[-1] for v in sp[1:])

                    # print sequence in FASTA format
                    print(">{}.{}.{}".format(i, j, k))
                    print_seq(first, subseq)

                except nx.NetworkXNoPath:
                    # there is no path, so do nothing
                    pass

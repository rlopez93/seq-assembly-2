#! /usr/bin/env python
from __future__ import print_function
from euler import find_eulerian_path
import sys
import screed
import networkx as nx

if len(sys.argv) != 3:
    print("Usage:", sys.argv[0], "<fast[aq]> <ksize>")
    sys.exit(1)

infilename = sys.argv[1]
k = int(sys.argv[2])

DG = nx.DiGraph()

for record in screed.open(infilename):
    seq = record['sequence']

    for i in xrange(len(seq)-k-1):
        kmer_a = seq[i:i+k]
        kmer_b = seq[i+1:i+k+1]
        DG.add_edge(kmer_a, kmer_b)


# function for printing sequence
# with 60 bases per row
def print_seq(a, seq):
    i = 0
    for base in a:
        print(base, sep='', end='')
        if (i+1) % 60 == 0:
            print()
        i += 1

    for base in seq:
        print(base, sep='', end='')
        if (i+1) % 60 == 0:
            print()
        i += 1
    print()

wc_subgraphs = nx.weakly_connected_component_subgraphs(DG)

for i, subgraph in enumerate(wc_subgraphs):
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
        # odd_cnt = 0
        in_deg = subgraph.in_degree
        out_deg = subgraph.out_degree
        sources = []
        targets = []

        for node in subgraph:
            ind = in_deg(node)
            outd = out_deg(node)
            # if ind != outd:
                # odd_cnt += 1
            if ind == 0:
                sources.append(node)
            if outd == 0:
                targets.append(node)

        paths = 0
        for j, source in enumerate(sources):
            for k, target in enumerate(targets):
                try:
                    sp = nx.shortest_path(subgraph, source, target)
                    first = sp[0]
                    subseq = (v[-1] for v in sp[1:])

                    print(">{}.{}.{}".format(i, j, k))
                    print_seq(first, subseq)
                except nx.NetworkXNoPath:
                    pass

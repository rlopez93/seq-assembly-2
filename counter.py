#! /usr/bin/env python
from __future__ import print_function
from collections import Counter
import screed
import sys

if len(sys.argv) != 4:
    print("Usage:", sys.argv[0], "<fast[aq]> <k> <C>")
    sys.exit(1)

infilename = sys.argv[1]
k = int(sys.argv[2])
C = int(sys.argv[3])

count = Counter()

for record in screed.open(infilename):
    seq = record['sequence']

    for i in xrange(len(seq)-k):
        count[seq[i:i+k]] += 1

for record in screed.open(infilename):
    seq = record['sequence']

    kmer = seq[:k]
    if len(seq) < k or count[kmer] < C:
        continue

    new_seq = list(kmer)
    i = 1
    while i+k <= len(seq):
        kmer = seq[i:i+k]
        if count[kmer] < C:
            break
        new_seq.append(kmer[-1])
        i += 1

    name = record['name']
    # annotations = record['annotations']
    accuracy = record['accuracy'][:len(new_seq)]

    print("@", name, sep='')
    print(''.join(new_seq))
    print('+\n{}'.format(accuracy))

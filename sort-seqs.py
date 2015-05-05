#! /usr/bin/env python
from __future__ import print_function
import screed
import sys

if len(sys.argv) != 2:
    print("Usage:", sys.argv[0], "<fast[aq]>")
    sys.exit(1)

seqs = []
for record in screed.open(sys.argv[1]):
    seqs.append(record['sequence'])


def print_seq(seq):
    for i, base in enumerate(seq):
        print(base, sep='', end='')
        if (i+1) % 60 == 0:
            print()
    print()

for i, seq in enumerate(sorted(seqs, key=len, reverse=True)):
    print(">", i, sep='')
    print_seq(seq)

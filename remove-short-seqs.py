#! /usr/bin/env python
from __future__ import print_function
import screed
import sys

if len(sys.argv) != 3:
    print("<py> <fast[aq]> <min_len>",
          file=sys.stderr)
    sys.exit(1)

infilename = sys.argv[1]
min_len = int(sys.argv[2])


def print_seq(seq):
    for i, base in enumerate(seq):
        print(base, end='')
        if (i+1) % 60 == 0:
            print()
    print()

for record in screed.open(infilename):
    name = record.name
    sequence = record.sequence

    if len(sequence) >= min_len:
        print(">", name, sep='')
        print_seq(sequence)

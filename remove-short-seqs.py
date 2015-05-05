#! /usr/bin/env python
import screed
import sys

if len(sys.argv) != 3:
    print >>sys.stderr, "<py> <fast[aq]> <min_len>"
    sys.exit(1)

infilename = sys.argv[1]
min_len = int(sys.argv[2])


def print_seq(seq):
    for i, base in enumerate(seq):
        sys.stdout.write(base)
        if (i+1) % 60 == 0:
            sys.stdout.write("\n")
    sys.stdout.write("\n")

for record in screed.open(infilename):
    name = record['name']
    sequence = record['sequence']

    if len(sequence) >= min_len:
        print ">{}".format(name)
        print_seq(sequence)

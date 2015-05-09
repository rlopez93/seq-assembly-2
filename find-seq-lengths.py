#! /usr/bin/env python2
from __future__ import print_function
import screed
import sys
from collections import Counter

if len(sys.argv) == 3:
    cutoff = int(sys.argv[2])
elif len(sys.argv) != 2:
    print("Usage:", sys.argv[0], "<fast[aq]>",
          file=sys.stderr)
    sys.exit(1)
else:
    cutoff = 500

infilename = sys.argv[1]

lengths = Counter()

for record in screed.open(infilename):
    curr_len = len(record['sequence'])
    curr_len -= curr_len % cutoff
    lengths[curr_len] += 1

print("   length  count")

for k, v in sorted(lengths.iteritems()):
    print("{0:4d}-{1:4d} {2:6d}".format(k, k+cutoff, v))


print("\nTotal contigs produced:", sum(lengths.values()))

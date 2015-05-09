#!/usr/bin/env python2
from __future__ import print_function
import fileinput

query = -1
length = 0
has_length = False
has_ident = False

good_contigs = []
bad_contigs = []

for line in fileinput.input():
    if line.startswith('Query='):
        query += 1
        has_ident = False
        has_length = False

    elif line.startswith('Length=') and not has_length:
        length = int(line.split('=')[1])
        has_length = True

    elif line.startswith(' Identities') and not has_ident:
        fields = line.split()
        matched = int(fields[2].split('/')[0])
        if matched != length:
            bad_contigs.append((query, length, matched))
        else:
            good_contigs.append((query, length, matched))
        has_ident = True

print("good contigs:", len(good_contigs), end='\n\n')
for query, length, matched in good_contigs:
    print(query, length, matched)
print()

print("bad contigs:", len(bad_contigs), end='\n\n')
for query, length, matched in bad_contigs:
    print(query, length, matched)
print()

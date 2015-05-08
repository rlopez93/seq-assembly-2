#! /usr/bin/env python
from __future__ import print_function
import fileinput

# coverage[i] -> coverage for gene at position i
coverage = []

for line in fileinput.input():
    fields = line.split()

    # the 8th column of fields contains the errors.
    # if fields has an 8th column, then we ignore it
    if len(fields) > 7:
        continue

    # get location of current read from map file
    location = int(fields[3])
    length = 36

    # calculate difference between end
    # of current read and size of coverage
    diff = location + length - len(coverage)

    # if current read extends past current coverage
    if diff > 0:
        # extend coverage to make up the difference
        coverage.extend([0] * diff)

    # update coverage
    for i in xrange(location, location+length):
        coverage[i] += 1

# now we find what intervals of the genome
# have non-zero coverage

i = 0
while i < len(coverage):
    # ignore positions with 0 coverage
    while i < len(coverage) and coverage[i] == 0:
        i += 1

    if i >= len(coverage):
        break

    # set start point for current coverage interval
    start = i
    sub_coverage = [coverage[i]]

    # find size of current coverage interval
    i += 1
    while i < len(coverage) and coverage[i] != 0:
        sub_coverage.append(coverage[i])
        i += 1

    total = len(sub_coverage)
    avg_cov = sum(sub_coverage) / float(total)

    print("{0:4d} {1:10d} {2:6.1f}"
          .format(total, start, avg_cov))

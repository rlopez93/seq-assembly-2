#! /usr/bin/env python
from __future__ import print_function
import fileinput

locations = []

for line in fileinput.input():
    fields = line.split()

    # # only error-free reads
    # if len(fields) > 7:
    #     continue

    location = int(fields[3])
    locations.append(location)

loc_max = max(locations)
loc_len = loc_max + 1 + 36

coverage = [0] * loc_len

for location in locations:
    for i in range(location, location+36):
        coverage[i] += 1

i = 0
while i < len(coverage):
    while i < len(coverage) and coverage[i] == 0:
        i += 1

    if i >= len(coverage):
        break

    start = i
    sub_coverage = [coverage[i]]

    i += 1
    while i < len(coverage) and coverage[i] != 0:
        sub_coverage.append(coverage[i])
        i += 1

    total = len(sub_coverage)
    avg_cov = sum(sub_coverage) / total

    print("{0:4d} {1:10d} {2:6.1f}"
          .format(total, start, avg_cov))

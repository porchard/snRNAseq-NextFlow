#!/usr/bin/env python
import os
import sys

infile = sys.argv[1]
counts = {}

with open(infile, 'r') as f:
    for line in f:
        line = line.rstrip()
        if line not in counts:
            counts[line] = 0
        counts[line] += 1

for line, count in counts.items():
    print('{}\t{}'.format(count, line))

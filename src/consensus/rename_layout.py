#!/usr/bin/env python

import sys
import re
from functools import partial

def read_map(fn):
    d = dict()
    with open(fn, "r") as f:
        for l in f:
            split = l.split()
            d['read'+split[0]] = ('read'+split[1])
    #return dict([(l.split()[1], l.split()[0]) for l in open(fn, "r")])
    return d

if len(sys.argv) < 2:
    print("Usage: %s <name map>.\nReads layout from standard input and prints converted layout to standard output." % sys.argv[0])
    sys.exit(1)

mapping = read_map(sys.argv[1])

for line in sys.stdin:
    fields = line.strip().split()
    print(fields[0], ','.join([(mapping[t[:-1]] + t[-1]) for t in fields[1].strip().split(',')]))

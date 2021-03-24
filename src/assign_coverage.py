#!/usr/bin/env python3

import sys
import re
from statistics import median_low
from collections import defaultdict

def parse_read_rec(r):
    return r.replace('read', '').replace('+', '').replace('-', '')

if len(sys.argv) < 3:
    print("Usage: %s <unitig layout> <read coverage>" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

layout_fn = sys.argv[1]
read_cov_fn = sys.argv[2]

read_cov=defaultdict(int)
with open(read_cov_fn, 'r') as f:
    for l in f:
        s = l.split()
        read_cov[s[0]] = int(s[1])

with open(layout_fn, 'r') as f:
    for l in f:
        s = l.split()
        print('%s %d' % (s[0], median_low([read_cov[parse_read_rec(r)] for r in s[1].split(',')])))

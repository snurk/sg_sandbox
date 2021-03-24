#!/usr/bin/env python3

import sys
import re
from statistics import median_low

def segment_length(split_line):
    if split_line[2] != '*':
        return len(split_line[2])
    else:
        for s in split_line:
            if s.startswith("LN:i:"):
                return int(s.split(":")[2])
        assert(False)

if len(sys.argv) < 3:
    print("Usage: %s <gfa> <segment coverage>" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

gfa_fn = sys.argv[1]
segment_cov_fn = sys.argv[2]

seg_cov=dict()
with open(segment_cov_fn, 'r') as f:
    for l in f:
        s = l.split()
        seg_cov[s[0]] = float(s[1])

with open(gfa_fn, 'r') as f:
    for l in f:
        if l.startswith('S\t'):
            s = l.split()
            seg = s[1]
            s = [item for item in s if not re.match("^(RC:i:|FC:i:|ll:f:).*", item)]
            length = segment_length(s)
            assert seg in seg_cov
            print('%s\tRC:i:%d' % ('\t'.join(s), int(length * seg_cov[seg])))
        else:
            print(l.strip())

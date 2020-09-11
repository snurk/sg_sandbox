#!/usr/bin/env python

import sys
import re
from statistics import median_low

if len(sys.argv) < 3:
    print("Usage: %s <gfa> <segment coverage>" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

gfa_fn = sys.argv[1]
segment_cov_fn = sys.argv[2]

seg_cov=dict()
with open(segment_cov_fn, 'r') as f:
    for l in f:
        s = l.split()
        seg_cov[s[0]] = int(s[1])

with open(gfa_fn, 'r') as f:
    for l in f:
        if l.startswith('S\t'):
            s = l.split()
            seg = s[1]
            assert s[3].startswith('LN')
            length = int(s[3].replace('LN:i:', ''))
            assert seg in seg_cov
            print('%s\tRC:i:%d' % (l.strip(), length * seg_cov[seg]))
        else:
            print(l.strip())

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

if len(sys.argv) < 3 or len(sys.argv) > 3 and sys.argv[3] != '--allow-absent':
    print("Usage: %s <gfa> <segment coverage> [--allow-absent]" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

allow_absent = False
gfa_fn = sys.argv[1]
segment_cov_fn = sys.argv[2]

if len(sys.argv) > 3:
    allow_absent = True

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
            assert seg in seg_cov or allow_absent
            cov = seg_cov[seg] if seg in seg_cov else 0.
            print('%s\tRC:i:%d' % ('\t'.join(s), int(length * cov)))
        else:
            print(l.strip())

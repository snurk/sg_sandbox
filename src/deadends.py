#!/usr/bin/env python3
import sys
#L       9806    -       9805    -

def segment_length(split_line):
    if split_line[2] != '*':
        return len(split_line[2])
    else:
        for s in split_line:
            if s.startswith("LN:i:"):
                return int(s.split(":")[2])
        assert(False)

def start_linked(o, first):
    return (o == '+') != first

if len(sys.argv) < 2:
    print("Usage: %s <in gfa> [length threshold = -1]" % sys.argv[0])
    print("Finds deadends")
    sys.exit(1)

len_thr = -1
if len(sys.argv) >= 3:
    len_thr = int(sys.argv[2])

with_start_l = set()
with_end_l = set()
with open(sys.argv[1], "r") as gfa:
    for l in gfa:
        if l.startswith("L"):
            split_col = l.split("\t")[1:5]
            if start_linked(split_col[1], True):
                with_start_l.add(split_col[0])
            else:
                with_end_l.add(split_col[0])

            if start_linked(split_col[3], False):
                with_start_l.add(split_col[2])
            else:
                with_end_l.add(split_col[2])

with open(sys.argv[1], "r") as gfa:
    for l in gfa:
        if l.startswith("S"):
            s = l.split("\t")
            length = segment_length(s)
            name = s[1]
            if len_thr > 0 and length > len_thr:
                continue
            if name not in with_start_l or name not in with_end_l:
                print(name)

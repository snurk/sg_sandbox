#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import sys
import re
#L       9806    -       9805    -

def length(split_line):
    if split_line[2] != '*':
        return len(split_line[2])
    else:
        for s in split_line:
            if s.startswith("LN:i:"):
                return int(s.split(":")[2])
        assert(False)

def read_lengths(gfa_fn):
    lengths = dict()
    with open(gfa_fn) as gfa:
        for l in gfa:
            if l[0] == 'S':
                sp=l.split('\t')
                lengths[sp[1]] = length(sp)
    return lengths

def ovl_len(cigar):
    match = re.match('^(\d+)M$', cigar)
    assert match
    return int(match.group(1))
    #pos = 0
    #for match in re.finditer('(\d+)M', cigar):
    #    pos += int(match.group(1))
    #for match in re.finditer('(\d+)I' if second else '(\d+)D', cigar):
    #    pos += int(match.group(1))
    #return pos

def inverse_sign(s):
    if s == "-":
        return "+"
    elif s == "+":
        return "-"
    else:
        assert(false)

def alt_enc(split):
    return "%s %s %s %s" % (split[2], inverse_sign(split[3]), split[0], inverse_sign(split[1]))

if len(sys.argv) < 2:
    print("Usage: %s <in gfa>" % sys.argv[0])
    print("Deduplicates copies of the same links. Outputs to stdout")
    sys.exit(1)

lengths = read_lengths(sys.argv[1])

used = dict()
with open(sys.argv[1], "r") as gfa:
    for l in gfa:
        if l.startswith("L"):
            split_col = l.split("\t")[1:6]
            link_str = " ".join(split_col[0:4])
            #try:
            ovl = ovl_len(split_col[4])
            if link_str not in used:
                ovl_bound = min(lengths[split_col[0]], lengths[split_col[2]]) - 1
                assert ovl_bound >= 0
                if ovl > ovl_bound:
                    print("For link '%s' overlap %d (cigar %s) will be trimmed to %d" % (link_str, ovl, split_col[4].rstrip(), ovl_bound), file=sys.stderr)
                    split_col[4] = ("%dM" % ovl_bound)
                    print('L\t%s' % '\t'.join(split_col))
                else:
                    print(l.rstrip())
            else:
                if used[link_str] != ovl:
                    print("For", link_str, "inconsistent overlaps:", used[link_str], "and", ovl, ". Diff", abs(ovl - used[link_str]), file=sys.stderr)
            used[alt_enc(split_col)] = ovl
            #bypassing a repeated-lines bug
            used[("\t".join(split_col))] = ovl
            #except ValueError:
            #   print("Oops, couldn't convert %s from line '%s'\n" % (split_col[4].strip()[:-1], l), file=sys.stderr)
            #   sys.exit(1)
        else:
            print(l.rstrip())

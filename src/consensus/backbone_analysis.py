#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from collections import defaultdict

import sys
import re
import numpy as np
import argparse

def swap(c):
    if c == '-':
        return '+'
    elif c == '+':
        return '-'
    assert False

def split_segment_desc(s_d):
    return (s_d[:-1], s_d[-1])

def transform_dir_seg(s):
    if s[0] == '>':
        return s[1:] + '+'
    elif s[0] == '<':
        return s[1:] + '-'
    assert(false)

parser = argparse.ArgumentParser(description="Join sequences of GFA nodes")
parser.add_argument("paths", help="Part of GAF detailing paths (name id1[+-],id2[+-],... start end)")
parser.add_argument("coverage", help="Node coverage file")
parser.add_argument("composition", help="Node composition stats")
#parser.add_argument("result", help="Output file")
args = parser.parse_args()

#if len(sys.argv) < 4:
#    print("Usage: script.py <contigs> <order_fn, with lines id1+,id2-,...> <gfa> <output fasta>")
#    exit(1)

#contigs_fn = sys.argv[1]
#order = sys.argv[2]
#out_fn = sys.argv[4]

seg_cov=dict()
with open(args.coverage, 'r') as f:
    for l in f:
        s = l.split()
        seg_cov[s[0]] = float(s[1])

seg_comp=dict()
with open(args.composition, 'r') as f:
    for l in f:
        s = l.split()
        seg_comp[s[0]] = int(s[1])

def match_length(a, i, j):
    assert a[i] == a[j]
    k = 1
    while i + k < len(a) and j + k < i:
        if a[i + k] != a[j + k]:
            break
        k += 1
    return k

def reuse_stats(segment_descs, excess_coeff, expected_cov):
    reused = list()
    bag = defaultdict(int)
    with_orientation = set()
    for s in segment_descs:
        assert s[:-1] in seg_cov
        assert s[:-1] in seg_comp
        if s[:-1] in bag:
            if not s in with_orientation:
                print("Trouble segment %s in the path in different orientations" % s[:-1], file=sys.stderr)
        #        assert False

        bag[s[:-1]] += 1
        with_orientation.add(s)

    for b in bag:
        if bag[b] > 1:
            print("Reuse of node %s (cov: %f, #reads: %d) %d times" % (b, seg_cov[b], seg_comp[b], bag[b]), file=sys.stderr)
            if seg_cov[b] < expected_cov * bag[b]:
                print("WARN", file=sys.stderr)
            if excess_coeff * bag[b] > seg_comp[b]:
                reused.append(b)
            else:
                print("No need to split at node %s: #reads %d, reused %d times" % (b, seg_comp[b], bag[b]), file=sys.stderr)

    return reused

def analyze_path(segment_descs):
    READ_EXCESS_COEFF = 2
    EXPECTED_COV = 15.
    reused = reuse_stats(segment_descs, READ_EXCESS_COEFF, EXPECTED_COV)
    return [i for i in range(0, len(segment_descs)) if segment_descs[i][:-1] in reused]
    #i = 0
    #n = len(segment_descs)
    #last_printed = 0
    #regions = list()
    #while i < n:
    #    s = segment_descs[i]
    #    j = 0
    #    minmax = 1000
    #    while j < i:
    #        s2 = segment_descs[j]
    #        if s[:-1] == s2[:-1]:
    #            if s[-1] != s2[-1]:
    #                print("Trouble segment %s in the path in different orientations (positions %d and %d)" % (s[:-1], i, j), file=sys.stderr)
    #                assert False
    #            else:
    #                print("Considering match %s between (positions %d and %d)" % (s, i, j), file=sys.stderr)
    #            m = match_length(segment_descs, i, j)
    #            print("Match length was %d" % m, file=sys.stderr)
    #            if minmax > m:
    #                print("Shorter than previous maximal match", file=sys.stderr)
    #                minmax = m
    #        j += 1

    #    if minmax != 1000:
    #        if i > last_printed:
    #            print("region %s" % ','.join(segment_descs[last_printed: i]), file=sys.stderr)
    #            regions.append(i)
    #        print("Will jump ahead %d positions" % minmax, file=sys.stderr)
    #        print("repetitive region %s" % ','.join(segment_descs[i: i + minmax]), file=sys.stderr)
    #        regions.append(i + minmax)
    #        last_printed = i + minmax

    #        i += minmax
    #    else:
    #        i += 1

    #if last_printed < n:
    #    print("region %s" % ','.join(segment_descs[last_printed: n]), file=sys.stderr)
    #    regions.append(n)
    #print("======================================", file=sys.stderr)

    #return regions

#with open(args.result, 'w') as out:
for line in open(args.paths, 'r'):
    print('=======================================', file=sys.stderr)
    s = line.split()
    l = s[1]
    #start = int(s[2])
    #end = int(s[3])
    #print('Processing path %s (start:end) -- %s (%d:%d)' % (s[0], l, start, end), file=sys.stderr)
    print('Processing path %s -- %s' % (s[0], l), file=sys.stderr)
    if len(l) == 0:
        continue

    if l[0] == '>' or l[0] == '<':
        segment_descs = [transform_dir_seg(t) for t in directed_seg_pattern.findall(l)]
    else:
        segment_descs = [t.strip() for t in l.split(',')]

    regions = analyze_path(segment_descs)
    regions.append(len(segment_descs))
    block = list()
    block_end = 0
    block_cnt = 1
    for boundary in regions:
        add = True
        for i in range(block_end, boundary):
            if segment_descs[i] in block:
                add = False
                break

        if add:
            block.extend(segment_descs[block_end: boundary])
        else:
            assert block_end > 0
            print('%s.%d %s' % (s[0], block_cnt, ','.join(block)))
            block_cnt += 1
            block = segment_descs[block_end: boundary]

        block_end = boundary

    if block_cnt == 1:
        print('%s %s' % (s[0], ','.join(block)))
    else:
        print('%s.%d %s' % (s[0], block_cnt, ','.join(block)))

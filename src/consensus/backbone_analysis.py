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

def swap_segment(s):
    return s[:-1] + swap(s[-1])

def split_segment_desc(s_d):
    return (s_d[:-1], s_d[-1])

def transform_dir_seg(s):
    if s[0] == '>':
        return s[1:] + '+'
    elif s[0] == '<':
        return s[1:] + '-'
    assert(false)

parser = argparse.ArgumentParser(description="Join sequences of GFA nodes")
parser.add_argument("paths", help="Part of GAF detailing paths (name id1[+-],id2[+-],...)")
parser.add_argument("coverage", help="Node coverage file")
parser.add_argument("composition", help="Node composition stats")
parser.add_argument("--excess-coeff", type=float, default=1., help="Do not split at nodes with >= (node-multiplicity * excess-coeff) reads available (default 1.)")
parser.add_argument("--min-cov", type=float, default=15., help="Minimal single-copy node coverage (default 15.)")
args = parser.parse_args()

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

def reuse_stats(segment_descs):
    reused = list()
    seg_mults = defaultdict(int)
    with_orientation = set()
    for s in segment_descs:
        s_n = s[:-1]
        assert s_n in seg_cov
        assert s_n in seg_comp
        if s_n in seg_mults:
            if not s in with_orientation:
                print("Trouble segment %s in the path in different orientations" % s_n, file=sys.stderr)
        #        assert False

        seg_mults[s_n] += 1
        with_orientation.add(s)

    for b in seg_mults:
        if seg_mults[b] > 1:
            if seg_cov[b] < args.min_cov * seg_mults[b]:
                print("WARN: Node %s (cov: %f, #reads: %d) used %d times" % (b, seg_cov[b], seg_comp[b], seg_mults[b]), file=sys.stderr)

            if seg_comp[b] < args.excess_coeff * seg_mults[b] - 0.01:
                print("Reuse of node %s (cov: %f, #reads: %d) %d times" % (b, seg_cov[b], seg_comp[b], seg_mults[b]), file=sys.stderr)
                reused.append(b)
            else:
                print("No need to split at node %s: #reads %d, reused %d times" % (b, seg_comp[b], seg_mults[b]), file=sys.stderr)

    return reused

def analyze_path(segment_descs):
    reused = reuse_stats(segment_descs)
    return [i for i in range(0, len(segment_descs)) if segment_descs[i][:-1] in reused]

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

    boundaries = analyze_path(segment_descs)
    boundaries.append(len(segment_descs))

    curr_block = list()
    curr_block_end = 0
    block_cnt = 1
    for boundary in boundaries:
        chunk = segment_descs[curr_block_end : boundary]
        print('Considering addition of chunk [%s] to the current block [%s]' % (','.join(chunk), ','.join(curr_block)), file=sys.stderr)
        #add = True
        #for i in range(curr_block_end, boundary):
        #    if segment_descs[i] in curr_block:
        #        add = False
        #        break
        if segment_descs[curr_block_end] not in curr_block and swap_segment(segment_descs[curr_block_end]) not in curr_block:
            print('Will extend', file=sys.stderr)
            curr_block.extend(chunk)
        else:
            print('Can not extend. Reporting current block as [%s.%d %s]' % (s[0], block_cnt, ','.join(curr_block)), file=sys.stderr)
            assert curr_block_end > 0
            print('%s.%d %s' % (s[0], block_cnt, ','.join(curr_block)))
            block_cnt += 1
            print('Resetting current block to', ','.join(chunk), file=sys.stderr)
            curr_block = chunk

        curr_block_end = boundary
        print('Updating block boundary to', curr_block_end, file=sys.stderr)

    if block_cnt == 1:
        print('%s %s' % (s[0], ','.join(curr_block)))
    else:
        print('%s.%d %s' % (s[0], block_cnt, ','.join(curr_block)))

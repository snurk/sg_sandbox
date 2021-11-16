#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import argparse

def swap(c):
    if c == '-':
        return '+'
    elif c == '+':
        return '-'
    assert False

def split_segment_desc(s_d):
    return (s_d[:-1], s_d[-1])

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
                l=re.sub('tig0+', '', l)
                sp=l.split('\t')
                lengths[sp[1]] = length(sp)
    return lengths

def read_cigars(gfa_fn):
    cigars = dict()
    with open(gfa_fn) as gfa:
        for l in gfa:
            if l[0] == 'L':
                l=re.sub('tig0+', '', l)
                sp=l.split('\t')
                cigars['%s%s|%s%s' % (sp[1], sp[2], sp[3], sp[4])] = sp[5]
                cigars['%s%s|%s%s' % (sp[3], swap(sp[4]), sp[1], swap(sp[2]))] = sp[5]
    return cigars

#def read_cigar(gfa_fn, s_d_1, s_d_2):
#    with open(gfa_fn) as gfa:
#        for l in gfa:
#            sp=l.split('\t')
#            if sp[0] == 'L' and re.sub('tig0+', '', sp[1]) == s_d_1[:-1] and re.sub('tig0+', '', sp[3]) == s_d_2[:-1]:
#                if sp[2] == s_d_1[-1] and sp[4] == s_d_2[-1]:
#                    return sp[5]
#
#            if sp[0] == 'L' and re.sub('tig0+', '', sp[3]) == s_d_1[:-1] and re.sub('tig0+', '', sp[1]) == s_d_2[:-1]:
#                if sp[4] == swap(s_d_1[-1]) and sp[2] == swap(s_d_2[-1]):
#                    return sp[5]
#    return None

def trim_len(cigar):#, second=True):
    match = re.match('^(\d+)M$', cigar)
    assert match
    return int(match.group(1))
    #pos = 0
    #for match in re.finditer('(\d+)M', cigar):
    #    pos += int(match.group(1))
    #for match in re.finditer('(\d+)I' if second else '(\d+)D', cigar):
    #    pos += int(match.group(1))
    #return pos

def transform_dir_seg(s):
    if s[0] == '>':
        return s[1:] + '+'
    elif s[0] == '<':
        return s[1:] + '-'
    assert(false)

def transform_back(s):
    if s[-1] == '+':
        return '>' + s[:-1]
    elif s[-1] == '-':
        return '<' + s[:-1]
    assert(false)

parser = argparse.ArgumentParser(description="Join sequences of GFA nodes")
parser.add_argument("paths", help="Part of GAF detailing paths (name [><]id1,[><]id2,... start end)")
parser.add_argument("gfa", help="GFA file specifying graph structure")
parser.add_argument("result", help="Output file")
parser.add_argument("--trusted-overhang", type=int, default=0, help="Required alignment past the first/last overlap")
parser.add_argument("--gaf-paths", action="store_true", help="Use GAF path format ([<>]-prefix instead of [+-]-suffix)")
args = parser.parse_args()

#if len(sys.argv) < 4:
#    print("Usage: script.py <contigs> <order_fn, with lines id1+,id2-,...> <gfa> <output fasta>")
#    exit(1)

#contigs_fn = sys.argv[1]
#order = sys.argv[2]
#out_fn = sys.argv[4]

print('Reading links from', args.gfa)
cigars = read_cigars(args.gfa)
lengths = read_lengths(args.gfa)

trusted_overhang_threshold = args.trusted_overhang

def get_cigar(s_d_1, s_d_2):
    q = s_d_1 + "|" + s_d_2
    return cigars[q] if q in cigars else None

def overlap_len(s_d_1, s_d_2):
    cigar = get_cigar(s_d_1, s_d_2)
    assert cigar is not None
    return trim_len(cigar)

def trim_start(segment_descs, start):
    if len(segment_descs) < 2:
        return 0, segment_descs

    ovl = overlap_len(segment_descs[0], segment_descs[1])
    assert start < lengths[segment_descs[0][:-1]]

    not_in_ovl = lengths[segment_descs[0][:-1]] - ovl
    if start >= not_in_ovl:
        print("Trimming redundant start of the path %s (path start: %d, overlap begins at %d)" % (segment_descs[0], start, not_in_ovl))
        assert not_in_ovl > 0
        return not_in_ovl, segment_descs[1:]
    elif not_in_ovl < trusted_overhang_threshold + start:
        assert not_in_ovl > 0
        print("Trimming untrusted start of the path %s (path start: %d, overlap begins at %d)" % (segment_descs[0], start, not_in_ovl))
        return not_in_ovl, segment_descs[1:]
    else:
        return 0, segment_descs

def total_len(segment_descs):
    if len(segment_descs) == 0:
        return 0

    assert segment_descs[0][:-1] in lengths
    tot = lengths[segment_descs[0][:-1]]

    for i in range(0, len(segment_descs)-1):
        ovl = overlap_len(segment_descs[i], segment_descs[i+1])
        print("Overlap between %s and %s is %d" % (segment_descs[i], segment_descs[i+1], ovl))
        assert segment_descs[i+1][:-1] in lengths
        tot += lengths[segment_descs[i+1][:-1]] - ovl

    return tot

def trim_end(segment_descs, tot, end):
    if len(segment_descs) < 2:
        return 0, segment_descs

    ovl = overlap_len(segment_descs[-2], segment_descs[-1])
    assert segment_descs[-1][:-1] in lengths
    assert end + lengths[segment_descs[-1][:-1]] > tot

    not_in_ovl = lengths[segment_descs[-1][:-1]] - ovl
    if tot - end >= not_in_ovl:
        print("Trimming redundant end of the path %s (path end: %d, path total: %d, last not in ovl %d)" % (segment_descs[-1], end, tot, not_in_ovl))
        return not_in_ovl, segment_descs[:-1]
    elif not_in_ovl < trusted_overhang_threshold + tot - end:
        print("Trimming untrusted end of the path %s (path end: %d, path total: %d, last not in ovl %d)" % (segment_descs[-1], end, tot, not_in_ovl))
        return not_in_ovl, segment_descs[:-1]
    else:
        return 0, segment_descs


directed_seg_pattern=re.compile(r"[<>][\w\-\_]+")
records = []

with open(args.result, 'w') as out:
    for line in open(args.paths, 'r'):
        print('=======================================')
        s = line.split()
        l = s[1]
        start = int(s[2])
        end = int(s[3])
        print('Processing path %s (start:end) -- %s (%d:%d)' % (s[0], l, start, end))
        if len(l) == 0:
            continue

        if l[0] == '>' or l[0] == '<':
            segment_descs = [transform_dir_seg(t) for t in directed_seg_pattern.findall(l)]
        else:
            segment_descs = [t.strip() for t in l.split(',')]

        while True:
            trimmed, segment_descs = trim_start(segment_descs, start)
            if trimmed == 0:
                break
            else:
                if start > trimmed:
                    start -= trimmed
                else:
                    print("Moving start")
                    start = 0

                assert end > trimmed
                end -= trimmed

        tot = total_len(segment_descs)
        #print("Total length", tot)
        while True:
            trimmed, segment_descs = trim_end(segment_descs, tot, end)
            if trimmed == 0:
                break
            else:
                assert tot > trimmed
                tot -= trimmed
                if end > tot:
                    print("Moving end")
                    end = tot

        if args.gaf_paths:
            print("%s %s %d %d" % (s[0], ''.join([transform_back(s) for s in segment_descs]), start, end), file=out)
        else:
            print("%s %s %d %d" % (s[0], ','.join(segment_descs), start, end), file=out)

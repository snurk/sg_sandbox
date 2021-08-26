#!/usr/bin/env python
import sys
import argparse
import parasail
import re
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def my_segment(s, flag):
    if flag == '-':
        return s.reverse_complement()
    elif flag == '+':
        return s
    else:
        assert False
        return ''

parser = argparse.ArgumentParser(description="Find overlaps and join contigs according to plan")
parser.add_argument("plan", help="File with join plan. Line format: 'name id1[+-],id2[+-],...'")
parser.add_argument("contigs", nargs='?', help="Input FASTA file (by default reading FASTA from stdin)")
parser.add_argument("-f", "--flank", type=int, default=25000, help="Flank region size to consider (default: 25000bp)")
parser.add_argument("-p", "--penalty", type=int, default=300, help="Mismatch/indel penalty (default: 300bp)")
parser.add_argument("-o", "--output", help="Output FASTA (by default prints to stdout)")
parser.add_argument("--min-ovl", type=int, default=0, help="Minimal length of overlap (default 0)")
parser.add_argument("--use-slow", action="store_true", help="Use slow DP implemenation")
parser.add_argument("--trim-both", action="store_true", help="Trim both 'first' and 'second' overlapping sequence (default: trim only second one)")
args = parser.parse_args()

if args.contigs:
    print("Reading from file", args.contigs, file=sys.stderr)
    i_handle = open(args.contigs, 'r')
else:
    print("Reading from stdin", file=sys.stderr)
    i_handle = sys.stdin

order = args.plan

if args.output:
    print("Writing to file", args.output, file=sys.stderr)
    o_handle = open(args.output, 'w')
else:
    print("Writing to stdin", file=sys.stderr)
    o_handle = sys.stdout

flank = args.flank
penalty = args.penalty

print('Will consider flanks of size', flank, file=sys.stderr)
print('Mismatch/indel penalty', penalty, file=sys.stderr)

if args.use_slow:
    print('Using slow DP implementation', file=sys.stderr)
else:
    print('Using Parasail', file=sys.stderr)

if args.trim_both:
    print('Will trim both joined sequences', file=sys.stderr)
else:
    print('Will trim only second joined sequence', file=sys.stderr)

def align_ends(s1, s2):
    ma = 1
    mm = -penalty
    ind = -penalty

    n=len(s1)
    m=len(s2)
    A = np.tile(np.arange(m + 1) * ind, (n+1, 1))

    for i in range(0, n+1):
        for j in range(0, m+1):
            m_s = mm
            if s1[i - 1] == s2[j - 1]:
                m_s = ma
            A[i, j] = max(A[i,j], A[i-1, j] + ind)
            A[i, j] = max(A[i,j], A[i, j-1] + ind)
            A[i, j] = max(A[i,j], A[i-1, j-1] + m_s)

    max_score=0
    s2_end_match=0
    for j in range(0, m+1):
        if A[n, j] > max_score:
            max_score = A[n, j]
            s2_end_match = j

    return s2_end_match, max_score

def parasail_ends(s1, s2):
    m = parasail.matrix_create("ACGT", 1, -penalty)
    res = parasail.sg_qb_de(s1, s2, penalty, penalty, m)
    return res.end_ref + 1, res.score

def overlap_and_join(s1, s2):
    f1 = min(len(s1), flank)
    f2 = min(len(s2), flank)
    #print('Aligning', s1[len(s1) - f1:], s2[:f2])
    print('Aligning', file=sys.stderr)
    if args.use_slow:
        pos, score = align_ends(s1[-f1:], s2[:f2])
    else:
        pos, score = parasail_ends(str(s1[-f1:]), str(s2[:f2]))

    assert score > 0
    assert pos > 4000
    print('Overlap', pos, file=sys.stderr)
    print('Score', score, file=sys.stderr)

    if args.trim_both:
        return s1[:-(pos // 2)] + s2[(pos // 2) + (pos % 2):]
    else:
        return s1 + s2[pos:]

def make_sequence(contig_dict, segment_order):
    segment_descs = [s.strip() for s in segment_order.split(',')]
    segments=[my_segment(contig_dict[s_d[:-1]], s_d[-1]) for s_d in segment_descs]

    print("Initiating with", segment_descs[0], file=sys.stderr)
    final = segments[0]
    for i in range(1, len(segments)):
        print("Appending", segment_descs[i], file=sys.stderr)
        final = overlap_and_join(final, segments[i])

    return final

contig_dict = dict()
for r in SeqIO.parse(i_handle, "fasta"):
    contig_dict[re.sub('tig0+', '', r.name)] = r.seq

if args.contigs:
    i_handle.close()

records = []
i = 0
for line in open(order, 'r'):
    i += 1
    print('=======================================', file=sys.stderr)
    s = line.split()
    l = s[1]
    print('Processing ', l, file=sys.stderr)
    if len(l) == 0:
        continue
    seq = make_sequence(contig_dict, l)
    name = s[0]
    print('Writing', name, file=sys.stderr)
    records.append(SeqRecord(seq, id=name, description=''))

SeqIO.write(records, o_handle, 'fasta')

if args.output:
    o_handle.close()

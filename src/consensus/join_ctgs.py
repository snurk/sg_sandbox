#!/usr/bin/env python
import sys
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
        sys.exit(1)
        return ''

def swap(c):
    if c == '-':
        return '+'
    elif c == '+':
        return '-'
    assert False

def split_segment_desc(s_d):
    return (s_d[:-1], s_d[-1])

def align_ends(s1, s2, ma = 1, mm = -300, ind = -300):
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

if len(sys.argv) < 4:
    print("Usage: script.py <contigs> <order_fn, with lines id1+,id2-,...> <output fasta> [flank, default: 25000bp][mm/ind penalty, default: 300]")
    exit(1)

contigs_fn = sys.argv[1]
order = sys.argv[2]
out_fn = sys.argv[3]

flank = 25000

if len(sys.argv) > 4:
    flank = int(sys.argv[4])

penalty = 300

if len(sys.argv) > 5:
    penalty = int(sys.argv[5])

print('Will consider flanks of size', flank)

def overlap_and_trim(s1, s2):
    f1 = min(len(s1), flank)
    f2 = min(len(s2), flank)
    #print('Aligning', s1[len(s1) - f1:], s2[:f2])
    print('Aligning')
    pos, score = align_ends(s1[len(s1) - f1:], s2[:f2], 1, mm = -penalty, ind = -penalty)
    print('Will trim', pos)
    print('Score', score)

    return s2[pos:]

def make_sequence(contig_dict, segment_order):
    segment_descs = [s.strip() for s in segment_order.split(',')]
    segments=[my_segment(contig_dict[s_d[:-1]], s_d[-1]) for s_d in segment_descs]

    final = segments[0]
    for i in range(0, len(segments)-1):
        print("Aligning", segment_descs[i], "and", segment_descs[i+1])
        final += overlap_and_trim(segments[i], segments[i+1])

    return final

contig_dict = dict()
for r in SeqIO.parse(open(contigs_fn, "r"), "fasta"):
    contig_dict[re.sub('tig0+', '', r.name)] = r.seq

records = []
i = 0
for line in open(order, 'r'):
    i += 1
    print('=======================================')
    s = line.split()
    l = s[1]
    print('Processing ', l)
    if len(l) == 0:
        continue
    seq = make_sequence(contig_dict, l)
    name = s[0]
    print('Writing', name)
    records.append(SeqRecord(seq, id=name, description=''))

SeqIO.write(records, open(out_fn, 'w'), 'fasta')

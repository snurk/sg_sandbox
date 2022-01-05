#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import argparse
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

def get_cigar(cigars, s_d_1, s_d_2):
    q = s_d_1 + "|" + s_d_2
    return cigars[q] if q in cigars else None

def trim_len(cigar, second=True):
    pos = 0
    for match in re.finditer('(\d+)M', cigar):
        pos += int(match.group(1))
    for match in re.finditer('(\d+)I' if second else '(\d+)D', cigar):
        pos += int(match.group(1))
    return pos

def make_sequence(cigars, contig_dict, segment_order, trim_flanks_to):
    segment_descs = [s.strip() for s in segment_order.split(',')]
    segments=[my_segment(contig_dict[s_d[:-1]], s_d[-1]) for s_d in segment_descs]

    tot = 0
    if trim_flanks_to < 0:
        final = segments[0]
        tot += len(segments[0])
    else:
        final = segments[0][-trim_flanks_to:]
        tot += trim_flanks_to
        print('Here taking', trim_flanks_to)

    for i in range(0, len(segments)-1):
        cigar = get_cigar(cigars, segment_descs[i], segment_descs[i+1])
        print(segment_descs[i], segment_descs[i+1])
        if (cigar is None):
           print("Cigar not known between %s and %s\n"%(segment_descs[i], segment_descs[i+1]))
        assert cigar is not None
        print('================== Basic step')
        print(segment_descs[i], segment_descs[i+1], cigar)
        pos = trim_len(cigar)
        print('From cigar', pos)
        l_before = tot
        final += segments[i+1][pos:]
        tot += (len(segments[i+1]) - pos)
        print('Added coordinates: [%d - %d)' % (l_before, tot))

    if trim_flanks_to > 0 and len(segments[-1]) > trim_flanks_to:
        final = final[:-(len(segments[-1]) - trim_flanks_to)]
        print('Here trimming', len(segments[-1]) - trim_flanks_to)

    return final, False
    #cigar = read_cigar(sys.argv[3], segment_descs[-1], segment_descs[0])
    #if cigar is not None:
    #    print('================== Circular step')
    #    print(segment_descs[-1], segment_descs[0], cigar)
    #    pos = trim_len(cigar)
    #    print('From cigar', pos)
    #    final = final[pos:]
    #    return final, True
    #else:
    #    return final, False

def transform_dir_seg(s):
    if s[0] == '>':
        return s[1:] + '+'
    elif s[0] == '<':
        return s[1:] + '-'
    assert(false)

parser = argparse.ArgumentParser(description="Join sequences of GFA nodes")
parser.add_argument("segments", help="FASTA file with graph segment sequences")
parser.add_argument("paths", help="File detailing paths (name id1[+-],id2[+-],...)")
parser.add_argument("gfa", help="GFA file specifying graph structure")
parser.add_argument("result", help="Output FASTA file")
parser.add_argument("--trim-flanks-to", type=int, default=-1, help="Trim flanking segments to specified value (disabled if negative)")
args = parser.parse_args()

if args.trim_flanks_to > 0:
    print('Will trim flanking segments to %dbp' % args.trim_flanks_to)

print('Collecting sequences from', args.segments)
contig_dict = dict()
for r in SeqIO.parse(open(args.segments, "r"), "fasta"):
    contig_dict[re.sub('tig0+', '', r.name)] = r.seq

print('Reading links from', args.gfa)
cigars = read_cigars(args.gfa)

directed_seg_pattern=re.compile(r"[<>]\w+")
records = []
for line in open(args.paths, 'r'):
    print('=======================================')
    s = line.split()
    if s[0] == 'name':
        continue

    l = s[1]
    print('Processing ', l)
    assert(len(l) != 0)

    if l[0] == '>' or l[0] == '<':
        l = ','.join([transform_dir_seg(s) for s in directed_seg_pattern.findall(l)])

    seq, circular = make_sequence(cigars, contig_dict, l, args.trim_flanks_to)
    name = s[0]
    print('Writing', name)
    records.append(SeqRecord(seq, id=name, description=''))

SeqIO.write(records, open(args.result, 'w'), 'fasta')

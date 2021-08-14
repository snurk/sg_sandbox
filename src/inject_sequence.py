#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Integrate node sequences into GFA")
parser.add_argument("segments", help="FASTA file with graph segment sequences")
parser.add_argument("gfa", nargs='?', help="GFA file specifying graph structure (by default reading from stdin)")
args = parser.parse_args()

print('Collecting sequences from', args.segments, file=sys.stderr)
contig_dict = dict()
for r in SeqIO.parse(open(args.segments, "r"), "fasta"):
    contig_dict[r.name] = r.seq

if args.gfa:
    print("Reading graph from", args.gfa, file=sys.stderr)
    stream = open(args.gfa, 'r')
else:
    print("Reading graph from stdin", file=sys.stderr)
    stream = sys.stdin

for l in stream:
    if l.startswith('S\t'):
        s = l.split()
        assert s[2] == '*' and s[1] in contig_dict
        s[2] = str(contig_dict[s[1]])
        print('\t'.join(s))
    else:
        print(l.strip())

if args.gfa:
    stream.close()

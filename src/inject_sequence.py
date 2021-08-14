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
parser.add_argument("gfa", help="GFA file specifying graph structure")
parser.add_argument("segments", help="FASTA file with graph segment sequences")
args = parser.parse_args()

print('Collecting sequences from', args.segments, file=sys.stderr)
contig_dict = dict()
for r in SeqIO.parse(open(args.segments, "r"), "fasta"):
    contig_dict[re.sub('tig0+', '', r.name)] = r.seq

with open(args.gfa, 'r') as f:
    for l in f:
        if l.startswith('S\t'):
            s = l.split()
            assert s[2] == '*' and s[1] in contig_dict
            s[2] = str(contig_dict[s[1]])
            print('\t'.join(s))
        else:
            print(l.strip())

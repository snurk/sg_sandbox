#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import os
from Bio import SeqIO

if len(sys.argv) < 3:
    print("Usage: %s <contigs_file> <output folder> [min_length = 0]" % sys.argv[0])
    sys.exit(1)

f_n = sys.argv[1]
folder = sys.argv[2]

min_len = 0

if len(sys.argv) > 3:
    min_len = int(sys.argv[3])

import os
if not os.path.exists(folder):
    os.makedirs(folder)

with open(f_n, "r") as ctgs_f:
    for ctg in SeqIO.parse(ctgs_f, "fasta"):
        if len(ctg.seq) > min_len:
            with open(os.path.join(folder, ctg.id + ".fasta"), "w") as out_f:
                SeqIO.write([ctg], out_f, "fasta")

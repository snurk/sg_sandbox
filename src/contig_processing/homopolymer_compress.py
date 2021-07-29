#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np

def compress_str(s):
    a = np.array(s)
    c = a[a != np.append('X', a[:-1])]
    return ''.join(c)

def compress_rec(r):
    r.seq = Seq(compress_str(r.seq))
    return r

compressed_seq_iterator = (compress_rec(record) for record in SeqIO.parse(sys.stdin, "fasta") if len(record.seq) > 0)
SeqIO.write(compressed_seq_iterator, sys.stdout, "fasta")

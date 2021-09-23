#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import select

def compress_str(s):
    a = np.array(s)
    c = a[a != np.append('X', a[:-1])]
    return ''.join(c)

def compress_rec(r):
    r.seq = Seq(compress_str(r.seq))
    return r

# Check that there is input ready
if sys.stdin in select.select([sys.stdin], [], [], 0.)[0]:
    compressed_seq_iterator = (compress_rec(record) for record in SeqIO.parse(sys.stdin, "fasta") if len(record.seq) > 0)
    SeqIO.write(compressed_seq_iterator, sys.stdout, "fasta")
else:
    print("Usage: %s reads fasta from stdin and writes to stdout" % sys.argv[0], file=sys.stderr)
    print("Tried reading from stdin, but it was empty!", file=sys.stderr)
    sys.exit(1)

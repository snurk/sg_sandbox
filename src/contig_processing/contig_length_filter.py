#!/usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO

if len(sys.argv) < 2:
    print("Usage: %s <length threshold> [FASTA file (default: stdin)]" % sys.argv[0])
    sys.exit(1)

i_fn = None
len_thr = int(sys.argv[1])
if len(sys.argv) >= 3:
    i_fn = sys.argv[2]

print("Retaining contigs of length >= ", len_thr, file=sys.stderr)

if i_fn:
    print("Reading from file", i_fn, file=sys.stderr)
    i_handle = open(i_fn, 'r')
else:
    print("Reading from stdin", file=sys.stderr)
    i_handle = sys.stdin

filtered_iterator = (record for record in SeqIO.parse(i_handle, "fasta") \
                      if len(record.seq) >= len_thr)

SeqIO.write(filtered_iterator, sys.stdout, "fasta")

if i_fn:
    i_handle.close()

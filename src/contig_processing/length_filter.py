#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import argparse
from Bio import SeqIO
import select

parser = argparse.ArgumentParser(description="Filter FASTA file based on length range and optionally trim sequences")
parser.add_argument("fasta", nargs='?', help="FASTA file (by default reading from stdin)")
parser.add_argument("--min", type=int, default=1, help="Minimal allowed length (default: 1)")
parser.add_argument("--max", type=int, default=sys.maxsize, help="Maximal allowed length (default: maxsize)")
parser.add_argument("--trim", action="store_true", help="Trim reads longer than maximal allowed length")
args = parser.parse_args()

print("Retaining contigs of length within [%d, %d]" % (args.min, args.max), file=sys.stderr)

if args.fasta:
    print("Reading from file", args.fasta, file=sys.stderr)
    i_handle = open(args.fasta, 'r')
else:
    # Check that there is input ready
    if sys.stdin in select.select([sys.stdin], [], [], 0.)[0]:
        print("Reading from stdin", file=sys.stderr)
        i_handle = sys.stdin
    else:
        parser.print_help(sys.stderr)
        print("Tried reading from stdin, but it was empty!", file=sys.stderr)
        sys.exit(1)

def check_len(n):
    return n >= args.min and (args.trim or n <= args.max)

def transform(r):
    if len(r.seq) > args.max:
        r.seq = r.seq[:args.max]
    return r

filtered_iterator = (transform(r) for r in SeqIO.parse(i_handle, "fasta") \
                      if check_len(len(r.seq)))

SeqIO.write(filtered_iterator, sys.stdout, "fasta")

if args.fasta:
    i_handle.close()

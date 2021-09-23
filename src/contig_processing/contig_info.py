#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import argparse
from Bio import SeqIO
import select

parser = argparse.ArgumentParser(description="Compute contig stats")
parser.add_argument("contigs", nargs='?', help="FASTA file with contigs (by default reading from stdin)")
parser.add_argument("--bed-mode", action="store_true", help="Generate valid bed file with records spanning entire contigs")
parser.add_argument("--coverage-keyword", help="Sequence used to mark the coverage in contig name, followed by a number (e.g. 'cov_')")
args = parser.parse_args()

def parse_coverage(ctg_name):
    m = re.search(args.coverage_keyword + "([\d\.]+)", ctg_name)
    if m:
        return float(m.group(1))
    else:
        return -1.

def count_gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

#print("Producing contig summary for", sys.argv[1])

if args.contigs:
    print("Reading from file", args.contigs, file=sys.stderr)
    i_handle = open(args.contigs, 'r')
else:
    # Check that there is input ready
    if sys.stdin in select.select([sys.stdin], [], [], 0.)[0]:
        print("Reading from stdin", file=sys.stderr)
        i_handle = sys.stdin
    else:
        parser.print_help(sys.stderr)
        print("Tried reading from stdin, but it was empty!", file=sys.stderr)
        sys.exit(1)

if not args.bed_mode:
    if args.coverage_keyword:
        print("name\tlength\tGC\tcoverage")
    else:
        print("name\tlength\tGC")

for record in SeqIO.parse(i_handle, "fasta"):
    name = record.name
    length = len(record.seq)
    if args.bed_mode:
        print("%s\t0\t%d" % (name, length))
    else:
        gc = count_gc(record.seq)
        if args.coverage_keyword:
            cov = parse_coverage(record.name)
            print("%s\t%d\t%f\t%f" % (name, length, gc, cov))
        else:
            print("%s\t%d\t%f" % (name, length, gc))

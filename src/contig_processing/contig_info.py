#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import gzip
from mimetypes import guess_type
from Bio import SeqIO

def gz_open(fn, mode):
    encoding = guess_type(fn)[1]  # uses file extension
    if encoding is None:
        print("Working with text file", fn)
        return open(fn, mode=mode)
    elif encoding == 'gzip':
        print("Working with gzipped file", fn)
        if mode == 'r':
            mode = 'rt'
        elif mode == 'w':
            mode = 'wt'
        else:
            raise ValueError('Unknown mode "{}"'.format(mode))
        return gzip.open(fn, mode=mode)
    else:
        raise ValueError('Unknown file encoding: "{}"'.format(encoding))

def parse_coverage(ctg_name):
    m = re.search("cov_([\d\.]+)", ctg_name)
    if m:
        return float(m.group(1))
    else:
        return -1.

def count_gc(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

if len(sys.argv) < 3:
    print("Usage: %s <contigs_file> <output>" % sys.argv[0])
    sys.exit(1)

#print("Producing contig summary for", sys.argv[1])

with open(sys.argv[2], "w") as output:
    output.write("name\tlength\tcoverage\tGC\n")
    with gz_open(sys.argv[1], 'r') as i_handle:
        for record in SeqIO.parse(i_handle, "fasta"):
            name = record.name
            length = len(record.seq)
            cov = parse_coverage(record.name)
            gc = count_gc(record.seq)
            output.write("%s\t%d\t%f\t%f\n" % (name, length, cov, gc))

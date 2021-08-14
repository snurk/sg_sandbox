#!/usr/bin/env python
from __future__ import print_function
import sys
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

if len(sys.argv) < 4:
    print("Usage: %s <length threshold> <contigs_file> <output>" % sys.argv[0])
    sys.exit(1)

len_thr = int(sys.argv[1])
i_fn = sys.argv[2]
o_fn = sys.argv[3]

print("Copying contigs of length >= ", len_thr, "from", i_fn, "to", o_fn)

with gz_open(i_fn, 'r') as i_handle:
    input_seq_iterator = SeqIO.parse(i_handle, "fasta")
    filtered_iterator = (record for record in input_seq_iterator \
                          if len(record.seq) >= len_thr)

    with gz_open(o_fn, 'w') as o_handle:
        SeqIO.write(filtered_iterator, o_handle, "fasta")

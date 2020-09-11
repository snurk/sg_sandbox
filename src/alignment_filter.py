#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import os
import sys
import re
import argparse
from pysam import AlignmentFile
import math

def guess_sam_bam(fn):
    if fn.endswith(".sam"):
        return ""
    elif fn.endswith(".bam"):
        return "b"
    else:
        raise ValueError("Unknown file extension")

def hard_clipped(cigartuples):
    prefix = 0
    suffix = 0
    assert cigartuples is not None
    #Hard-clip -- code 5
    if len(cigartuples) > 0 and cigartuples[0][0] == 5:
        prefix = cigartuples[0][1]
    if len(cigartuples) > 1 and cigartuples[-1][0] == 5:
        suffix = cigartuples[-1][1]
    return prefix, suffix

def insertion_cnt(cigartuples):
    total = 0
    assert cigartuples is not None
    #Insertion -- code 1
    for code, cnt in cigartuples:
        if code == 1:
            total += cnt
    return total

parser = argparse.ArgumentParser(description="Filter alignments based on various criteria")
parser.add_argument("sam",
        help="Input SAM/BAM alignment file")
parser.add_argument("--min-len", type=int, default=0, help="Minimal length of alignment (default 0)")
parser.add_argument("--query-frac", type=float, default=0., help="Minimal fraction of the query sequence covered (default 0.)")
parser.add_argument("--min-idy", type=float, default=0., help="Minimal alignment identity defined as 1. - (edit_distance / alignment_matrix_length) (default 0.)")
parser.add_argument("--use-secondary", action="store_true", help="Minimal alignment identity")
parser.add_argument("--use-supplementary", action="store_true")
parser.add_argument("--use-all", action="store_true")
parser.add_argument("--filtered", help="Output SAM/BAM file")
parser.add_argument("--bed", help="Output bed file")
parser.add_argument("--bed-trim", type=int, default=0, help="Cut this many bp from each side of the covered region")
parser.add_argument("--overhang-margin", type=int, default=-1, help="Margin for over-hang reads, negative to disable over-hangs (default: -1)")
args = parser.parse_args()

print("Filtering alignments from", args.sam, file=sys.stderr)
alignment = AlignmentFile(args.sam, 'r' + guess_sam_bam(args.sam), check_sq=False)

if args.query_frac > 1.:
    print("Minimal covered fraction of the query sequence can not exceed 1. Specified", args.query_frac, file=sys.stderr)
    sys.exit(1)

if args.min_idy > 1.:
    print("Minimal alignment identity can not exceed 1. Specified", args.min_idy, file=sys.stderr)
    sys.exit(2)

print("Using min alignment length %d, min query fraction %.2f and min alignment identity %.2f" % (args.min_len, args.query_frac, args.min_idy), file=sys.stderr)

if args.bed:
    print("Will write covered region info to BED file", args.bed, file=sys.stderr)
    print("Will trim ", args.bed_trim, " bp from each side of covered region")
    out_bed = open(args.bed, 'w')

if args.filtered:
    print("Will write filtered alignments to", args.filtered, file=sys.stderr)
    out_alignment = AlignmentFile(args.filtered, 'w' + guess_sam_bam(args.filtered), template=alignment)

#FIXME improve logic
if args.overhang_margin >= 0:
    print("Will allow overhangs with tolerance of ", args.overhang_margin, "bp", file=sys.stderr)

min_len = args.min_len
overhang_margin = args.overhang_margin
query_frac = args.query_frac
min_idy = args.min_idy

def identity(align_len, edit_dist):
    return 1. - (edit_dist / align_len)

def check_query_cov(r, query_length, ref_length):
    assert query_length > 0
    if (r.query_alignment_length / query_length) > query_frac:
        return True

    if overhang_margin >= 0:
        if r.query_alignment_end + overhang_margin > query_length and r.reference_start < overhang_margin:
            return True
        if r.query_alignment_start < overhang_margin and r.reference_end + overhang_margin > ref_length:
            return True

    return False

print("qname\tqlen\tqstart\tqend\treverse\trname\trlen\trstart\trend")

for r in alignment:
    if r.is_unmapped:
        continue

    if r.is_secondary and not args.use_secondary and not args.use_all:
        continue

    if r.is_supplementary and not args.use_supplementary and not args.use_all:
        continue

    assert r.cigarstring is not None

    if r.query_alignment_length < min_len:
        continue

    if min_idy > 0.:
        assert r.has_tag('NM')

    if r.has_tag('NM'):
        #idy=identity(r.reference_length, r.get_tag('NM'))
        idy=identity(r.reference_length + insertion_cnt(r.cigartuples), r.get_tag('NM'))
        #idy=identity(r.query_alignment_length, r.get_tag('NM'))
    else:
        idy = 1.

    if idy < min_idy:
        continue

    init_length = r.infer_read_length()
    if not check_query_cov(r, init_length, alignment.lengths[r.reference_id]):
        continue

    if r.is_reverse:
        clipped_suffix, _ = hard_clipped(r.cigartuples)
        query_alignment_start = init_length - r.query_alignment_end - clipped_suffix
        query_alignment_end = init_length - r.query_alignment_start - clipped_suffix
        reference_start_str = "(%d" % r.reference_end
        reference_end_str = "%d]" % r.reference_start
    else:
        clipped_prefix, _ = hard_clipped(r.cigartuples)
        query_alignment_start = r.query_alignment_start + clipped_prefix
        query_alignment_end = r.query_alignment_end + clipped_prefix
        reference_start_str = "[%d" % r.reference_start
        reference_end_str = "%d)" % r.reference_end

    print("Alignment of query %s (length %d) : [%d - %d) to %s (length %d) : %s - %s. Identity: %.3f. Length on query (reference): %d (%d). Reverse: %r" % \
            (r.query_name, init_length, query_alignment_start, query_alignment_end, \
                alignment.get_reference_name(r.reference_id), alignment.lengths[r.reference_id], \
                reference_start_str, reference_end_str, \
                idy * 100., r.query_alignment_length, r.reference_length, r.is_reverse), \
            file=sys.stderr)

    print("%s\t%d\t%d\t%d\t%r\t%s\t%d\t%d\t%d" % (r.query_name, init_length, query_alignment_start, query_alignment_end, r.is_reverse, alignment.get_reference_name(r.reference_id), alignment.lengths[r.reference_id], r.reference_start, r.reference_end))

    if args.filtered:
        out_alignment.write(r)

    if args.bed:
        if r.reference_length < 2 * args.bed_trim:
            print("WARN reference span %d to small for used trim %d" % (r.reference_length, args.bed_trim), file=sys.stderr)
        else:
            print("%s\t%d\t%d\t%s" % (alignment.get_reference_name(r.reference_id), r.reference_start + args.bed_trim, r.reference_end - args.bed_trim, r.query_name), file=out_bed)

if args.filtered:
    out_alignment.close()

if args.bed:
    out_bed.close()

alignment.close()

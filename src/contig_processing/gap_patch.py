#!/usr/bin/env python

import os
import sys
import re
import argparse
import math
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(description="Create (raw) file for gap patching")
parser.add_argument("ctg_info", help="Contig length info file")
parser.add_argument("align_info", help="Alignment info file")
parser.add_argument("--max-trim", type=int, default=200, help="Maximal number of base-pairs allowed to be trimmed from contigs (default: 200)")
parser.add_argument("--reliable-len", type=int, default=10000, help="Minimal alignment length to be considered (default: 10000)")
parser.add_argument("--max-gap", type=int, default=10000, help="Maximal gap length that can be patched (default: 10000)")
args = parser.parse_args()

print("Reading contig details from", args.ctg_info, file=sys.stderr)

ctgs = pd.read_table(args.ctg_info)

print("Reading alignment details from", args.align_info, file=sys.stderr)

df = pd.read_table(args.align_info)

df.sort_values(by=['qname', 'qstart'])

#df = df[(df.qstart < 100) | ((df.qlen - df.qend) < 100)]
#print(df.head())

#def get_start_patch(row):
#    if not row.reverse:
#        return "%s\t%d\t%d\t.\t0\t+" % (row.rname, 0, row.rstart), row.rstart
#    else:
#        return "%s\t%d\t%d\t.\t0\t-" % (row.rname, row.rend, row.rlen), row.rlen - row.rend
#
#def get_end_patch(row):
#    if not row.reverse:
#        return "%s\t%d\t%d\t.\t0\t+" % (row.rname, row.rend, row.rlen), row.rlen - row.rend
#    else:
#        return "%s\t%d\t%d\t.\t0\t-" % (row.rname, 0, row.rstart), row.rstart

contig_start_mappings = defaultdict(list)
contig_end_mappings = defaultdict(list)

MAX_TRIM = args.max_trim
RELIABLE_LEN = args.reliable_len
MAX_GAP = args.max_gap

start_map = dict()
end_map = dict()

for _, row in ctgs.iterrows():
    #row.name returns number of the row, so using row[0] instead
    start_map[row[0]] = 0
    end_map[row[0]] = row.length

def report_ctg(name):
    assert start_map[name] < end_map[name]
    print("%s ctg: %s\t%d\t%d\t.\t0\t+" % (name, name, start_map[name], end_map[name]))

def compatible(e, s):
    assert e.rname == s.rname
    closing = False
    if e.reverse != s.reverse:
        return False
    if not e.reverse:
        if s.rstart > e.rend and s.rstart < e.rend + MAX_GAP:
            print("Closing gap between %s and %s with %s+" % (e.qname, s.qname, e.rname), file=sys.stderr)
            print("%s patch: %s\t%d\t%d\t.\t0\t+" % (e.qname, e.rname, e.rend, s.rstart))
            return True
    else:
        if e.rstart > s.rend and e.rstart < s.rend + MAX_GAP:
            print("Closing gap between %s and %s with %s-" % (e.qname, s.qname, e.rname), file=sys.stderr)
            print("%s patch: %s\t%d\t%d\t.\t0\t-" % (e.qname, s.rname, s.rend, e.rstart))
            return True
    return False

for _, row in df.iterrows():
    #if len(curr_qname) > 0 and not curr_reported and start_patched:
    #    report_ctg(curr_qname)

    #curr_qname = row.qname
    #start_patched = False
    #end_patched = False
    #curr_reported = False
    #start_map[row.qname] = 0
    #end_map[row.qname] = row.qlen

    if row.rlen < RELIABLE_LEN:
        continue

    if row.qstart < MAX_TRIM:
        contig_start_mappings[row.rname].append(row)

    if row.qlen - row.qend < MAX_TRIM:
        contig_end_mappings[row.rname].append(row)

for ref in contig_start_mappings:
    if ref not in contig_end_mappings:
        continue

    print("Considering ref tig", ref, file=sys.stderr)
    for e in contig_end_mappings[ref]:
        for s in contig_start_mappings[ref]:
            if compatible(e, s):
                if e.qend != e.qlen:
                    print("Setting nontrivial ending point for %s, trimming %d" % (e.qname, e.qlen - e.qend), file=sys.stderr)
                if s.qstart != 0:
                    print("Setting nontrivial starting point for %s, trimming %d" % (s.qname, s.qstart), file=sys.stderr)
                end_map[e.qname] = e.qend
                start_map[s.qname] = s.qstart

if len(start_map) == 1:
    print("Nothing to process", file=sys.stderr)

for c in start_map:
    report_ctg(c)

#if not curr_reported and start_patched:
#    report_ctg(curr_qname)

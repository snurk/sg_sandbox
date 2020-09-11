#!/usr/bin/env python

import os
import sys
import re
import argparse
import math
import pandas as pd

parser = argparse.ArgumentParser(description="Create bed file for telomere patching")
parser.add_argument("align_info", help="Input alignment file")
args = parser.parse_args()

print("Reading alignment details from", args.align_info, file=sys.stderr)

df = pd.read_table(args.align_info)

df.sort_values(by=['qname', 'qstart'])

#df = df[(df.qstart < 100) | ((df.qlen - df.qend) < 100)]
#print(df.head())

start_map = dict()
end_map = dict()

def report_ctg(name):
    print("%s\t%d\t%d\t.\t0\t+" % (name, start_map[name], end_map[name]))

def get_start_patch(row):
    if not row.reverse:
        return "%s\t%d\t%d\t.\t0\t+" % (row.rname, 0, row.rstart), row.rstart
    else:
        return "%s\t%d\t%d\t.\t0\t-" % (row.rname, row.rend, row.rlen), row.rlen - row.rend

def get_end_patch(row):
    if not row.reverse:
        return "%s\t%d\t%d\t.\t0\t+" % (row.rname, row.rend, row.rlen), row.rlen - row.rend
    else:
        return "%s\t%d\t%d\t.\t0\t-" % (row.rname, 0, row.rstart), row.rstart

curr_qname = ""
start_patched = False
end_patched = False
curr_reported = False

for index, row in df.iterrows():
    if row.qname != curr_qname:
        assert row.qname not in start_map
        assert row.qname not in end_map

        if len(curr_qname) > 0 and not curr_reported and start_patched:
            report_ctg(curr_qname)

        curr_qname = row.qname
        start_patched = False
        end_patched = False
        curr_reported = False
        start_map[row.qname] = 0
        end_map[row.qname] = row.qlen
    else:
        assert row.qname in start_map
        assert row.qname in end_map

    if row.qstart < 100:
        patch, patch_len = get_start_patch(row)
        if patch_len > row.qstart + 100:
            assert not start_patched
            start_map[row.qname] = row.qstart
            start_patched = True
            print("Successfully patching start of %s by %d bp" % (row.qname, patch_len), file=sys.stderr)
            print(patch)

    if row.qlen - row.qend < 100:
        patch, patch_len = get_end_patch(row)
        if patch_len > row.qlen - row.qend + 100:
            assert not end_patched
            end_map[row.qname] = row.qend
            end_patched = True
            assert not curr_reported
            report_ctg(curr_qname)
            curr_reported = True
            print("Successfully patching end of %s by %d bp" % (row.qname, patch_len), file=sys.stderr)
            print(patch)

if not curr_reported and start_patched:
    report_ctg(curr_qname)

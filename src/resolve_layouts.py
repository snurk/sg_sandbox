#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import re
import sys
from collections import defaultdict
import argparse

gaf_path_pattern=re.compile(r"[<>][\w\-\_]+")

def transform_from_gaf(s):
    if s[0] == '>':
        return s[1:] + '+'
    elif s[0] == '<':
        return s[1:] + '-'
    assert(false)

def transform_to_gaf(s):
    if s[-1] == '+':
        return '>' + s[:-1]
    elif s[-1] == '-':
        return '<' + s[:-1]
    assert(false)

def process_path(l):
    if l[0] == '>' or l[0] == '<':
        return [transform_from_gaf(t) for t in gaf_path_pattern.findall(l)]
    else:
        return [t.strip() for t in l.split(',')]

parser = argparse.ArgumentParser(description="Recursive resolution of the layouts")
parser.add_argument("fragment_names", help="File with names of the fragments to resolve OR layout file OR .gfa file -- segment names will be extracted")
parser.add_argument("composition", help="File specifying composition of intermediate fragments")
parser.add_argument("--resolved-marker", help="String marking the nodes that got duplicated during repeat resolution (e.g. '_i'). Part of the name past the marker will be ignored.")
parser.add_argument("--partial-resolve", action="store_true", help="Allow 'partial' resolving of layout. By default fails if can't resolve down to path of <resolved-prefix>.* fragments")
parser.add_argument("--gaf-paths", action="store_true", help="Use GAF path format ([<>]-prefix instead of [+-]-suffix)")
parser.add_argument("--miniasm", help="File with miniasm layout via 'a' lines")
parser.add_argument("--resolved-prefix", default='read', help="Prefix of completely resolved elements")
args = parser.parse_args()

mapping=defaultdict(list)

usage=defaultdict(int)

if args.miniasm:
    print("Loading miniasm layout", file=sys.stderr)
    with open(args.miniasm, 'r') as miniasm_layout:
        for l in miniasm_layout:
            if not l.startswith('a'):
                continue

            s = l.split()
            assert s[0] == 'a'

            mapping[s[1]].append(s[3] + s[4])
    print("Loaded", file=sys.stderr)

print("Loading compression mappings", file=sys.stderr)
with open(args.composition, 'r') as compress_mapping:
    for l in compress_mapping:
        s = l.split()
        assert s[0] not in mapping
        mapping[s[0]] = process_path(s[1])

print("Loaded", file=sys.stderr)

names = []

if args.fragment_names.endswith('.gfa'):
    print("Loading fragment names from GFA file", file=sys.stderr)
    with open(args.fragment_names, 'r') as gfa_fn:
        for l in gfa_fn:
            if not l.startswith('S'):
                continue

            #TODO optimize for GFA with sequences
            s = l.split()
            assert s[0] == 'S'
            names.append(s[1])
    print("Loaded", file=sys.stderr)
else:
    print("Loading fragment info from text file", file=sys.stderr)
    with open(args.fragment_names, 'r') as fragment_handle:
        for l in fragment_handle:
            s = l.strip().split()
            names.append(s[0].strip())
            if len(s) > 1:
                mapping[s[0]] = process_path(s[1])

    print("Loaded", file=sys.stderr)

def opposite_sign(o):
    if o == '+':
        return '-'
    elif o == '-':
        return '+'
    assert False

def swap_sign(n_o):
    return n_o[:-1] + opposite_sign(n_o[-1])

def need_swap(n_o):
    if n_o[-1] == '+':
        return False
    elif n_o[-1] == '-':
        return True
    assert False

def resolve(n_o, resolved):
    if n_o.startswith(args.resolved_prefix):
        resolved.append(n_o)
        return

    name = n_o[:-1]
    if args.resolved_marker:
        pos=name.find(args.resolved_marker)
        assert pos != 0
        if pos > 0:
            name=name[:pos]

    if not args.partial_resolve and name not in mapping:
        print("Problem with ", name, file=sys.stderr)

    assert args.partial_resolve or name in mapping

    if usage[name] > 0:
        print(name, "used multiple times", file=sys.stderr)

    usage[name] += 1

    if name not in mapping and args.partial_resolve:
        resolved.append(name + n_o[-1])
        return

    parts = mapping[name]

    if need_swap(n_o):
        for i in range(len(parts),0,-1):
            p = parts[i - 1]
            resolve(swap_sign(p), resolved)
    else:
        for p in parts:
            resolve(p, resolved)

print("Resolving segment composition", file=sys.stderr)
for s in names:
    resolved = []
    resolve(s + '+', resolved)
    if args.gaf_paths:
        print("%s\t%s" % (s, ''.join([transform_to_gaf(s) for s in resolved])))
    else:
        print("%s\t%s" % (s, ','.join(resolved)))

print("Layout resolution done", file=sys.stderr)

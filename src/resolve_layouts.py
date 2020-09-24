#!/usr/bin/env python3

import sys
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Recursive resolution of the layouts")
parser.add_argument("fragment_names", help="File with names of the fragments to resolve (or .gfa file -- segment names will be extracted)")
parser.add_argument("composition", help="File specifying composition of intermediate fragments")
parser.add_argument("--duplicate-shift", type=int, default=0, help="Ids shift to do if dealing with multiple occurences of the same fragment")
parser.add_argument("--partial-resolve", action="store_true", help="Allow 'partial' resolving of layout. By default fails if can't resolve down to path of readXXX fragments")
parser.add_argument("--miniasm", help="File with miniasm layout via 'a' lines")
args = parser.parse_args()

duplicate_shift = args.duplicate_shift

assert duplicate_shift == 0 or not args.partial_resolve

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
        mapping[s[0]] = s[1].split(',')
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
    print("Loading fragment names from text file", file=sys.stderr)
    with open(args.fragment_names, 'r') as names_handle:
        for l in names_handle:
            names.append(l.strip())
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

def resolve(n_o, resolved, duplicate = False):
    if n_o.startswith('read'):
        if not duplicate or duplicate_shift == 0:
            #minor optimization
            resolved.append(n_o)
        else:
            i = int(n_o[4:-1])
            print("Shifting %d into %d" % (i, i + duplicate_shift), file=sys.stderr)
            resolved.append('read%d%s' % ((i + duplicate_shift), n_o[-1]))
        return

    name = n_o[:-1]
    assert args.partial_resolve or name in mapping
    if duplicate_shift > 0:
        assert usage[name] < 2

    assert not duplicate or usage[name] > 0
    d = usage[name] > 0
    #if d and duplicate_shift > 0:
    if d:
        print("Will use duplicates to resolve", name, file=sys.stderr)
        #assert False

    usage[name] += 1

    if name not in mapping and args.partial_resolve:
        resolved.append(n_o)
        return

    parts = mapping[n_o[:-1]]

    if need_swap(n_o):
        for i in range(len(parts),0,-1):
            p = parts[i - 1]
            resolve(swap_sign(p), resolved, d)
    else:
        for p in parts:
            resolve(p, resolved, d)

print("Resolving segment composition", file=sys.stderr)
for s in names:
    resolved = []
    resolve(s + '+', resolved)
    print(s, ','.join(resolved))

print("All done", file=sys.stderr)

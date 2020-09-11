#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) < 4:
    print("Usage: %s <in GFA> <miniasm GFA with layout ('a' records)> <compress mapping file>" % sys.argv[0], file=sys.stderr)
    print("Output to stdout, brief log to stderr")
    sys.exit(1)

#a       utg000001l      0       read47654       +       637

mapping = defaultdict(list)

print("Loading miniasm layout", file=sys.stderr)
with open(sys.argv[2], 'r') as miniasm_layout:
    for l in miniasm_layout:
        if not l.startswith('a'):
            continue

        s = l.split()
        assert s[0] == 'a'

        mapping[s[1]].append(s[3] + s[4])
print("Loaded", file=sys.stderr)

print("Loading compression mappings", file=sys.stderr)
with open(sys.argv[3], 'r') as compress_mapping:
    for l in compress_mapping:
        s = l.split()
        assert s[0] not in mapping
        mapping[s[0]] = s[1].split(',')
print("Loaded", file=sys.stderr)

segments = []
print("Loading graph segments", file=sys.stderr)
with open(sys.argv[1], 'r') as gfa_fn:
    for l in gfa_fn:
        if not l.startswith('S'):
            continue

        s = l.split()
        assert s[0] == 'S'
        segments.append(s[1])
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
    if n_o.startswith('read'):
        resolved.append(n_o)
        return

    assert n_o[:-1] in mapping

    parts = mapping[n_o[:-1]]

    if need_swap(n_o):
        for i in range(len(parts),0,-1):
            p = parts[i - 1]
            resolve(swap_sign(p), resolved)
    else:
        for p in parts:
            resolve(p, resolved)

print("Resolving segment composition", file=sys.stderr)
for s in segments:
    resolved = []
    resolve(s + '+', resolved)
    print(s, ','.join(resolved))

print("All done", file=sys.stderr)

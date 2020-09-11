#!/usr/bin/env python

import sys
import re

if len(sys.argv) < 2:
    print("Usage: %s <paf>" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

def modify(a, b, l, o, line):
    assert a == 0 or b == 0 or a == l or b == l, "here " + line
    if a == 0:
        return 0, o
    if b == 0:
        return o, 0
    if a == l:
        return l, l - o
    if b == l:
        return l - o, l

fn = sys.argv[1]
with open(fn, 'r') as f:
    for l in f:
        sl = l.split()
        assert(len(sl) == 13)
        #if sl[9] != sl[10]:
        #    print("help", sl[9], sl[10])
        #assert(sl[9] == sl[10])

        l1 = int(sl[1])
        a1 = int(sl[2])
        b1 = int(sl[3])

        l2 = int(sl[6])
        a2 = int(sl[7])
        b2 = int(sl[8])

        ovl1 = b1 - a1
        ovl2 = b2 - a2

        if ovl1 >= l1 or ovl2 >= l2:
            #FIXME check with Bri!!!
            print("Ovl size problem", l, file=sys.stderr)
            continue

        assert ovl1 < l1 and ovl2 < l2

        #FIXME check with Bri!!!
        if a1 != 0 and b1 != 0 and a1 != l1 and b1 != l1:
            print("Ignored 1 not extremal", l, file=sys.stderr)
            continue

        if a2 != 0 and b2 != 0 and a2 != l2 and b2 != l2:
            print("Ignored 2 not extremal", l, file=sys.stderr)
            continue

        ovl = min(ovl1, ovl2)

        a1, b1 = modify(a1, b1, l1, ovl, l)
        a2, b2 = modify(a2, b2, l2, ovl, l)

        sl[2] = str(a1)
        sl[3] = str(b1)

        sl[7] = str(a2)
        sl[8] = str(b2)

        print('\t'.join(sl))

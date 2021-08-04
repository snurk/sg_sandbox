#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
import sys

def reshuffle_element(in1, in2, pos):
    global global_cnt
    global reshuffle_flag
    global no_break
    assert len(in1) == len(in2)
    assert pos <= len(in1)
    if pos == len(in1):
        global_cnt += 1
        print("Alt%d.1" % global_cnt, ','.join(in1).replace(no_break, ','))
        print("Alt%d.2" % global_cnt, ','.join(in2).replace(no_break, ','))
        return
    if in1[pos] != in2[pos]:
        #prevent swapping of first differing position
        if not reshuffle_flag:
            reshuffle_flag = True
        else:
            in1_c = in1[:]
            in2_c = in2[:]
            in1_c[pos], in2_c[pos] = in2[pos], in1[pos]
            reshuffle_element(in1_c, in2_c, pos + 1)
    reshuffle_element(in1, in2, pos + 1)

global_cnt=0
reshuffle_flag=False

if len(sys.argv) < 3:
    print("Usage: %s <s1> <s2> (comma separated lists of the same size) [non-break symbol (default: ':')]" % sys.argv[0])
    sys.exit(1)

s1 = sys.argv[1]
s2 = sys.argv[2]

no_break = ':'

if len(sys.argv) > 3:
    no_break = sys.argv[3]

reshuffle_element(s1.split(','), s2.split(','), 0)

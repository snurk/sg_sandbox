#!/usr/bin/env python

import sys
import re


class UnitigFollowFailure(Exception):
    """Raised when couldn't follow unitig"""
    pass

def opposite_sign(o):
    if o == '+':
        return '-'
    elif o == '-':
        return '+'
    assert False

def swap_sign(n_o):
    return n_o[:-1] + opposite_sign(n_o[-1])

def swap_all(c):
    return [swap_sign(r) for r in reversed(c)]

if len(sys.argv) < 3:
    print("Usage: %s <original read backbone layout> <unitig layout>" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

backbone_layout_fn = sys.argv[1]
unitig_layout_fn = sys.argv[2]

read_2_unitig = dict()
unitig_2_reads = dict()

print("Reading unitig layout file", unitig_layout_fn, file=sys.stderr)
with open(unitig_layout_fn, 'r') as f:
    for l in f:
        name, layout_str = l.split()
        layout = [r for r in layout_str.split(',')]
        unitig_2_reads[name] = layout
        for pos, r in enumerate(layout):
            r_name = r[:-1]
            if r_name in read_2_unitig:
                assert False
                #assert read_2_unitig[r_name][0] == name
                #read_2_unitig[r_name][1] = -1
            read_2_unitig[r_name] = (name, pos)

def containing_unitig(r):
    r_name = r[:-1]
    if r_name not in read_2_unitig:
        return ""
    name, pos = read_2_unitig[r_name]
    r_in_layout = unitig_2_reads[name][pos]
    assert r_in_layout.startswith(r_name)
    if r_in_layout[-1] == r[-1]:
        return name + '+'
    else:
        return name + '-'

def following_unitig_layout(r):
    r_name = r[:-1]
    assert(r_name in read_2_unitig)
    name, pos = read_2_unitig[r_name]
    assert(pos >= 0)
    layout = unitig_2_reads[name]
    r_in_layout = layout[pos]
    assert r_in_layout.startswith(r_name)
    if r_in_layout[-1] == r[-1]:
        return layout[pos:]
    else:
        return swap_all(layout[:pos+1])

def analyze_concordance(contig_layout, i, unitig_layout):
    j = 0
    n = len(contig_layout)
    m = len(unitig_layout)

    assert i < n
    assert m > 0
    assert contig_layout[i] == unitig_layout[j]

    MAX_SKIP = 2

    #TODO use Smith-Waterman?
    while i < n and j < m:
        #print("Trying to skip reads absent in unitigs")
        while i < n and (contig_layout[i][:-1] not in read_2_unitig):
            print("Skipping contig-read %s absent in unitigs" % contig_layout[i][:-1], file=sys.stderr)
            i+=1

        if i == n:
            print("WARNING: Wasn't able to locate 'suffix' contig-read(s) in unitigs", file=sys.stderr)
            break

        skipped = 0
        while j < m and contig_layout[i] != unitig_layout[j]:
            print("Trying to skip unitig-read %s since it didn't match the contig" % unitig_layout[j], file=sys.stderr)
            j+=1
            skipped += 1

        if j == m:
            print("WARNING: Last read of the unitig was skipped", file=sys.stderr)
            #print("HERE: Current 'contig' read is ", contig_layout[i], file=sys.stderr)
            #return -1
            return -1 if skipped > MAX_SKIP else i

        assert contig_layout[i] == unitig_layout[j]

        i+=1
        j+=1

    #FIXME stub for more involved procedure

    #Came to the end of the current unitig
    if j == m:
        return i

    #Came to the end of the contig, unitig goes forward
    if i == n:
        print("Unitig goes past the end of the contig", file=sys.stderr)
        print("Extended by %d reads" % (m - j), file=sys.stderr)
        return i

    print("Stopped at j=%d (m=%d)" %(j, m), file=sys.stderr)
    print("PROBLEM: %s in contig did not match %s in unitig" % (contig_layout[i], unitig_layout[j]), file=sys.stderr)
    return -1

def examine_layout(contig_layout):
    layout_str = ""
    layout_delim = ""
    i = 0
    n = len(contig_layout)
    #flag showing that no read has been matched to the contig yet
    unmatched_prefix = True
    unmatched_suffix = False
    all_good = True
    while i < n:
        r = contig_layout[i]
        utg = containing_unitig(r)

        if utg == "":
            if unmatched_prefix:
                print("WARNING: Wasn't able to locate 'prefix' contig-read %s in unitigs" % (r), file=sys.stderr)
            else:
                print("WARNING: Wasn't able to locate contig-read %s in unitigs considering it as start of 'suffix'" % (r), file=sys.stderr)
                unmatched_suffix = True
            i += 1
        else:
            print("Next read %s was located within utg %s" % (r, utg), file=sys.stderr)
            try:
                if unmatched_suffix:
                    #reset flag
                    unmatched_suffix = False
                    #Something suddenly matched again
                    #print("Couldn't match some reads at the unitigs boundary")
                    raise UnitigFollowFailure("Couldn't match some reads at the unitigs boundary", file=sys.stderr)

                further_unitig_layout = following_unitig_layout(r)
                utg_start = (len(further_unitig_layout) == len(unitig_2_reads[utg[:-1]]))
                if not utg_start:
                    if unmatched_prefix:
                        print("Unitig %s layout extends to the left of the contig" % utg, file=sys.stderr)
                        print("Extended by %d reads" % (len(unitig_2_reads[utg[:-1]]) - len(further_unitig_layout)), file=sys.stderr)
                    else:
                        print("PROBLEM: Discordant starting points", file=sys.stderr)
                        #raise UnitigFollowFailure("Discordant starting points")

                unmatched_prefix = False
                i_jump = analyze_concordance(contig_layout, i, further_unitig_layout)
                if i_jump < 0:
                    #print("Discordant alignment, couldn't recover")
                    unmatched_prefix = True
                    assert containing_unitig(contig_layout[i]) == utg
                    print("PROBLEM: Skipping beyond problematig unitig ", utg, file=sys.stderr)
                    while i < n and (containing_unitig(contig_layout[i]) == "" or containing_unitig(contig_layout[i]) == utg):
                        i += 1
                    raise UnitigFollowFailure("Discordant alignment")
                else:
                    i = i_jump
            except UnitigFollowFailure as utg_fail:
                print("PROBLEM", utg_fail, file=sys.stderr)
                all_good = False
                layout_delim += "!!!"

            layout_str += (layout_delim + utg)
            layout_delim = ","

    #if i == n:
    #    return True
    #return False
    assert(i == n)
    return layout_str, all_good

with open(backbone_layout_fn, 'r') as f:
    for l in f:
        contig_name, contig_layout_str = l.split()
        contig_layout = contig_layout_str.split(',')
        print("Examining layout of contig,", contig_name, file=sys.stderr)
        layout_str, all_good = examine_layout(contig_layout)
        if all_good:
            print("Successfully uplifted the layout for contig", contig_name, file=sys.stderr)
        else:
            print("Couldn't uplift the layout for contig", contig_name, file=sys.stderr)
        print("%s %s" % (contig_name, layout_str))

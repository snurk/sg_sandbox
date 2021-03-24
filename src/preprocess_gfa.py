#!/usr/bin/env python3
import sys
#L       9806    -       9805    -

def inverse_sign(s):
    if s == "-":
        return "+"
    elif s == "+":
        return "-"
    else:
        assert(false)

def alt_enc(split):
    return "%s %s %s %s" % (split[2], inverse_sign(split[3]), split[0], inverse_sign(split[1]))

if len(sys.argv) < 2:
    print("Usage: %s <in gfa>" % sys.argv[0])
    print("Deduplicates copies of the same links. Outputs to stdout")
    sys.exit(1)

used = dict()
with open(sys.argv[1], "r") as gfa:
    for l in gfa:
        if l.startswith("L"):
            split_col = l.split("\t")[1:6]
            link_str = " ".join(split_col[0:4])
            ovl = int(split_col[4][:-1])
            if link_str not in used:
                print(l.rstrip())
            else:
                if used[link_str] != ovl:
                    print("For", link_str, "inconsistent overlaps:", used[link_str], "and", ovl, ". Diff", abs(ovl - used[link_str]), file=sys.stderr)
            used[alt_enc(split_col)] = ovl
            #bypassing a repeated-lines bug
            used[("\t".join(split_col))] = ovl
        else:
            print(l.rstrip())

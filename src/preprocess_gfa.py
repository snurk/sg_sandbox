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
    return "%s\t%s\t%s\t%s" % (split[2], inverse_sign(split[3]), split[0], inverse_sign(split[1]))

if len(sys.argv) < 2:
    print("Usage: %s <in gfa>" % sys.argv[0])
    print("Deduplicates copies of the same links. Outputs to stdout")
    sys.exit(1)

used = set()
with open(sys.argv[1], "r") as gfa:
    for l in gfa:
        if l.startswith("L"):
            split_col = l.split("\t")[1:5]
            if "\t".join(split_col) not in used:
                print(l.rstrip())
            used.add(alt_enc(split_col))
            #bypassing a repeated-lines bug
            used.add("\t".join(split_col))
        else:
            print(l.rstrip())

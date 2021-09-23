#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "script.sh <reads.fasta> <ovl.paf> [min_overlap = 2000]"
    exit 239
fi

scripts_root=$(dirname $(readlink -e $0))/..
reads=$1
ovl_paf=$2

echo "Will use reads from $reads"
echo "Will use overlaps from $ovl_paf"

if [ "$#" -lt 3 ]; then
    min_ovl=2000
else
    min_ovl=$3
    echo "Will use min_ovl threshold of $min_ovl"
fi

echo "Run microasm to build string graph"
#Run microasm to build string graph
#"Fixed" gfakludge is so bad that it can not process gfa files without header correctly
echo -e "H\tVN:Z:1.0" > microasm.gfa
$scripts_root/../miniasm/miniasm -e 0 -d 0 -n 0 -m $min_ovl -s $min_ovl -f $reads $ovl_paf >> microasm.gfa 2> microasm.log

# Deduplicate copies of same links.
# + Trick to make gfapy process the custom tag
echo "Deduplicating links"
$scripts_root/preprocess_gfa.py microasm.gfa 2> preprocess.log | sed 's/SD/sd/g' > processed.gfa

echo "Microasm construction complete"

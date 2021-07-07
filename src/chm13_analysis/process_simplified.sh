#!/bin/bash
set -eou

if [ "$#" -lt 4 ]; then
    echo "Usage: $(basename $0) <simplified.gfa> <init.gfa> <node_composition.txt> <read_coverage.txt> [compact_prefix=m_]"
    exit 1
fi

prefix="m_"
if [ "$#" -gt 4 ]; then
    prefix=$5
fi

scripts_root=$(dirname $(readlink -e $0))/..

simplified=$1
name=$(basename $simplified .gfa)
init=$2
resolved_mapping=$3
read_cov=$4

grep "^S" $simplified | cut -f 2 > $name.tmp
echo -e "H\tVN:Z:1.0" > $name.fixed.gfa
$scripts_root/../tangle-resolution/scripts/extract_subgraph.py $name.tmp < $init | grep "^S" | sed $'s/\tRC.*//g' >> $name.fixed.gfa
#grep -f $name.tmp $init | sed $'s/\tRC.*//g' >> $name.fixed.gfa
grep "^L" $simplified >> $name.fixed.gfa
#rm -f $name.tmp

$scripts_root/compact_gfa.py $name.fixed.gfa $name.tmp.gfa $prefix 2> $name.final_compress_mapping.txt

cat $resolved_mapping $name.final_compress_mapping.txt > $name.mapping.txt
$scripts_root/resolve_layouts.py $name.tmp.gfa $name.mapping.txt > $name.resolved_mapping.txt

$scripts_root/assign_coverage.py $name.resolved_mapping.txt $read_cov > $name.recompressed.cov
$scripts_root/inject_coverage.py $name.tmp.gfa $name.recompressed.cov > $name.recompressed.gfa

echo -e "H\tVN:Z:1.0" > $name.recompressed.noseq.gfa
$scripts_root/../gfacpp/gfatools/gfatools view -S $name.recompressed.gfa >> $name.recompressed.noseq.gfa

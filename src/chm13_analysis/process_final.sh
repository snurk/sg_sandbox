#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "Usage: $(basename $0) <graph.gfa> <node_composition.txt> <out_prefix> [compact_prefix=m_] [read_coverage.txt]"
    exit 1
fi

comp_prefix="m_"
if [ "$#" -gt 3 ]; then
    comp_prefix=$4
fi

scripts_root=$(dirname $(readlink -e $0))/..

graph=$1
resolved_mapping=$2
prefix=$3

mkdir -p $(dirname $prefix)

$scripts_root/compact_gfa.py <(sed 's/\tCL.*//g' $graph) $prefix.tmp.gfa $comp_prefix 2> $prefix.final_compress_mapping.txt

cat $resolved_mapping $prefix.final_compress_mapping.txt > $prefix.mapping.txt
#Trick to ignore bait nodes
grep -v f_ $prefix.tmp.gfa > $prefix.tmp2.gfa
$scripts_root/resolve_layouts.py $prefix.tmp2.gfa $prefix.mapping.txt --resolved-marker _i > $prefix.resolved_mapping.txt
#$scripts_root/resolve_layouts.py $prefix.tmp.gfa $prefix.mapping.txt --resolved-marker _i > $prefix.resolved_mapping.txt

if [ "$#" -gt 4 ]; then
    read_cov=$5
    $scripts_root/assign_coverage.py $prefix.resolved_mapping.txt $read_cov > $prefix.cov
    $scripts_root/inject_coverage.py $prefix.tmp.gfa $prefix.cov --allow-absent > $prefix.gfa
else
    cp $prefix.tmp.gfa $prefix.gfa
fi

echo -e "H\tVN:Z:1.0" > $prefix.noseq.gfa
$scripts_root/../gfacpp/gfatools/gfatools view -S $prefix.gfa >> $prefix.noseq.gfa

rm -f $prefix.tmp.gfa
rm -f $prefix.tmp2.gfa

echo "Processing done"

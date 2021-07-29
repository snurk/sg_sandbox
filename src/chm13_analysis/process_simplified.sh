#!/bin/bash
set -eou

if [ "$#" -lt 5 ]; then
    echo "Usage: $(basename $0) <graph.gfa> <init.gfa> <node_composition.txt> <read_coverage.txt> <out_prefix> [compact_prefix=m_]"
    exit 1
fi

comp_prefix="m_"
if [ "$#" -gt 5 ]; then
    comp_prefix=$6
fi

scripts_root=$(dirname $(readlink -e $0))/..

graph=$1
init=$2
resolved_mapping=$3
read_cov=$4
prefix=$5

mkdir -p $(dirname $prefix)

$scripts_root/fix_node_info.sh $graph $init $prefix.fixed

$scripts_root/compact_gfa.py $prefix.fixed.gfa $prefix.tmp.gfa $comp_prefix 2> $prefix.final_compress_mapping.txt

cat $resolved_mapping $prefix.final_compress_mapping.txt > $prefix.mapping.txt
$scripts_root/resolve_layouts.py $prefix.tmp.gfa $prefix.mapping.txt --resolved-marker _i > $prefix.resolved_mapping.txt

$scripts_root/assign_coverage.py $prefix.resolved_mapping.txt $read_cov > $prefix.recompressed.cov
$scripts_root/inject_coverage.py $prefix.tmp.gfa $prefix.recompressed.cov > $prefix.recompressed.gfa

echo -e "H\tVN:Z:1.0" > $prefix.recompressed.noseq.gfa
$scripts_root/../gfacpp/gfatools/gfatools view -S $prefix.recompressed.gfa >> $prefix.recompressed.noseq.gfa

rm $prefix.tmp
rm $prefix.tmp.gfa

echo "Processing done"

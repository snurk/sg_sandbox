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

echo "Processing resolved and/or manually pruned graph"

mkdir -p $(dirname $prefix)
rm -f $prefix.compact_mapping.txt

echo "Normalizing graph"
$scripts_root/preprocess_gfa.py $graph > $prefix.tmp.gfa

echo "Compacting graph"
$scripts_root/../gfacpp/build/test $prefix.tmp.gfa $prefix.tmp.gfa --compact --prefix $comp_prefix --id-mapping $prefix.compact_mapping.txt &> $prefix.compact.log
#$scripts_root/compact_gfa.py <(sed 's/\tCL.*//g' $graph) $prefix.tmp.gfa $comp_prefix 2> $prefix.compact_mapping2.txt > $prefix.compact2.log

#Trick to ignore bait nodes
#grep -v f_ $prefix.tmp.gfa > $prefix.tmp2.gfa
#$scripts_root/resolve_layouts.py $prefix.tmp2.gfa $prefix.mapping.txt --resolved-marker _i > $prefix.resolved_mapping.txt

echo "Inferring unitig backbones"
$scripts_root/resolve_layouts.py $prefix.tmp.gfa <(cat $resolved_mapping $prefix.compact_mapping.txt) --resolved-marker _i > $prefix.resolved_mapping.txt

if [ "$#" -gt 4 ]; then
    read_cov=$5
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py $prefix.resolved_mapping.txt $read_cov > $prefix.cov
    #Trick to ignore bait nodes
    #$scripts_root/inject_coverage.py --allow-absent $prefix.cov $prefix.tmp.gfa > $prefix.gfa
    $scripts_root/inject_coverage.py $prefix.cov $prefix.tmp.gfa > $prefix.gfa
else
    mv $prefix.tmp.gfa $prefix.gfa
fi

echo "Getting graph without sequences"
echo -e "H\tVN:Z:1.0" > $prefix.noseq.gfa
$scripts_root/../gfacpp/gfatools/gfatools view -S $prefix.gfa >> $prefix.noseq.gfa

rm -f $prefix.tmp*

echo "Processing done"

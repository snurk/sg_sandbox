#!/bin/bash
set -e

scripts_root=$(dirname $(readlink -e $0))/..

#processed=$1
#init=$2
processed=<...>
name=$(basename $processed .gfa)
init=<...>
resolved_mapping=../../resolved_mapping.txt
read_cov=../../min_read.cov

grep "^S" $processed | awk '{print $2}' | sed 's/^/S\\s/g' | sed 's/$/\\s/g' > $name.tmp
head -1 $init > $name.fixed.gfa
grep -f $name.tmp $init >> $name.fixed.gfa
#grep -f $name.tmp $init | sed $'s/\tRC.*//g' >> $name.fixed.gfa
grep "^L" $processed >> $name.fixed.gfa
#rm -f $name.tmp

conda activate
$scripts_root/compact_gfa.py $name.fixed.gfa $name.tmp.gfa 2> $name.final_compress_mapping.txt

cat $resolved_mapping $name.final_compress_mapping.txt > $name.mapping.txt
$scripts_root/resolve_layouts.py $name.tmp.gfa $name.mapping.txt > $name.resolved_mapping.txt

$scripts_root/assign_coverage.py $name.resolved_mapping.txt $read_cov > $name.recompressed.cov
$scripts_root/inject_coverage.py $name.tmp.gfa $name.recompressed.cov > $name.recompressed.gfa

echo -e "H\tVN:Z:1.0" > $name.recompressed.noseq.gfa
$scripts_root/../gfacpp/gfatools/gfatools view -S $name.recompressed.gfa >> $name.recompressed.noseq.gfa

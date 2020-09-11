#!/bin/bash
set -e

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
~/git/ngs_scripts/bogart_gfa/compress.py $name.fixed.gfa $name.tmp.gfa 2> $name.final_compress_mapping.txt

cat $resolved_mapping $name.final_compress_mapping.txt > $name.mapping.txt
touch stub.microasm.gfa
~/git/ngs_scripts/bogart_gfa/resolve_mapping.py $name.tmp.gfa stub.microasm.gfa $name.mapping.txt > $name.resolved_mapping.txt

~/git/ngs_scripts/gfakluge/assign_coverage.py $name.resolved_mapping.txt $read_cov > $name.recompressed.cov
~/git/ngs_scripts/gfakluge/inject_coverage.py $name.tmp.gfa $name.recompressed.cov > $name.recompressed.gfa

echo -e "H\tVN:Z:1.0" > $name.recompressed.noseq.gfa
~/git/gfatools/gfatools view -S $name.recompressed.gfa >> $name.recompressed.noseq.gfa

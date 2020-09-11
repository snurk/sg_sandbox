#!/bin/bash
set -e

#CHR=$1

## assume mapList is your current list of nodes
CHR=chr16

grep "^S" ../../8k/chr_analysis/$CHR.*.noseq.gfa | awk '{print $2}' | sort | uniq > mapList

cat ../output.filtered.vcf | grep BND | grep -v chr |awk '{split($NF, a, ":"); if (a[3] > 1) print $0}' > $CHR.vcf

java -classpath /data/korens/devel/utils/ SubFile mapList $CHR.vcf |awk '{print $1" "$5}' |tr 'N' ' ' | tr '[' ' ' |tr ']' ' ' | tr ':' ' '|awk '{print $1; print $2}'|sort |uniq > list

echo "Sanity check"
java -classpath /data/korens/devel/utils/ SubFile list ../../8k/simplified.nodes.hg38.out |awk '{if ($4-$3 > 50000 && $NF > 99) print $6}'|sort |uniq -c

cat list mapList |sort |uniq > tmp
mv tmp list

diff list mapList

#cp /data/Phillippy/t2t-share/assemblies/drafts/20200602/simplified.noseq.gfa ./
#/home/nurks2/git/ngs_scripts/gfakluge/neighborhood simplified.noseq.gfa $CHR.gfa list 10000

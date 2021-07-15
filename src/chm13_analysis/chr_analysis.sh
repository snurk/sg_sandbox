#!/bin/bash
set -e

module load bedtools

#~/git/gfatools/gfatools view -S simplified.gfa > simplified.noseq.gfa
root=$(dirname $(readlink -e $0))

out_dir=$1
mashmap=$(readlink -e $2)
gfa=$(readlink -e $3)
gfa_noseq=$(readlink -e $4)
colors=$(readlink -e $5)
chrinfo=$(readline -e $6)

mkdir -p $out_dir
cd $out_dir

ln -sf $colors
ln -sf $chrinfo

# && $10 > 95.
#Filter mappings
awk '{if ($4 > 500000 + $3) print $0}' $mashmap | sort -k1,1 > filtered.out

#Find nodes associated with a single chromosome
awk '{print $1,$6}' filtered.out | sort | uniq | awk '{print $1}' | sort -k1,1 | uniq -u > good.nodes.txt

#Take info about only good nodes
join good.nodes.txt filtered.out > good.nodes.out

awk '{print $6,$8,$9,$1,$10}' good.nodes.out | sed 's/ /\t/g' > good.nodes.bed
bedtools coverage -a $chrinfo -b good.nodes.bed | sort -k 7 -n -r > frac.txt

#under interactive job
#~/useful_scripts/dump_coverage.sh assembly 0.00001 min_cov.txt

#~/git/ngs_scripts/gfakluge/assign_coverage.py ../resolved_mapping.txt ../min_read.cov > simplified.cov

rm -f chr*.txt chr*.log chr*.gfa
for chr in $(awk {'print $1'} $colors) ; do
    echo "Processing $chr"
    frac=$(grep "${chr}\s" frac.txt | awk '{print $7}' | grep -Po "\.\\d\\d" | sed 's/\.//g')
    grep "${chr}\s" good.nodes.out | awk '{print $1}' | sort | uniq > $chr.txt
    #$root/../../gfacpp/build/neighborhood $gfa $chr.$frac.gfa $chr.txt 10000 &> $chr.log
    $root/../../gfacpp/build/neighborhood $gfa_noseq $chr.$frac.noseq.gfa $chr.txt 10000 &> $chr.noseq.log
done

echo "Name,color,chr" > color.csv
#get color.csv
join <(awk '{print $6,$1}' good.nodes.out | sort | uniq | sort -k1,1) <(sort -k1,1 $colors) | awk '{print $2,$3,$1}' | sed 's/ /,/g' >> color.csv

cd -

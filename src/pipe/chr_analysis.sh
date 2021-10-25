#!/bin/bash
set -eou

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <out_dir> <mashmap.out> <gfa> <chr.colors> <chr.bed> [alignment size threshold, default=50000] [identity threshold, default=98]"
    exit 1
fi

root=$(dirname $(readlink -e $0))

out_dir=$1
mashmap=$(readlink -e $2)
gfa=$(readlink -e $3)
colors=$(readlink -e $4)
chrinfo=$(readlink -e $5)

aln_size=50000
if [ "$#" -gt 5 ]; then
    aln_size=$6
fi

identity=98
if [ "$#" -gt 6 ]; then
    identity=$7
fi

mkdir -p $out_dir
cd $out_dir

# && $10 > 95.
#Filter mappings
awk '{if ($4 >= '$aln_size' + $3 && $NF >= '$identity') print $0}' $mashmap | sort -k1,1 > filtered.out

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
$root/../extract_cov_micro.sh < $gfa > coverage.csv

for chr in $(awk {'print $1'} $colors) ; do
    echo "Processing $chr"
    frac=$(grep "${chr}\s" frac.txt | awk '{print $7}' | grep -Po "\.\\d\\d" | sed 's/\.//g')
    grep "${chr}\s" good.nodes.out | awk '{print $1}' | sort | uniq > $chr.txt
    $root/../../gfacpp/build/neighborhood $gfa $chr.$frac.gfa -n $chr.txt -r 10000 -c coverage.csv &> $chr.neib.log
done

echo "Name,color,chr" > color.csv
#get color.csv
join <(awk '{print $6,$1}' good.nodes.out | sort | uniq | sort -k1,1) <(sort -k1,1 $colors) | awk '{print $2,$3,$1}' | sed 's/ /,/g' >> color.csv

cd -

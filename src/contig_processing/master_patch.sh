#!/bin/bash
set -eou

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <chr_name> <reference> <fragments> <out_dir> [threads, default = 8]"
    exit 1
fi

root=$(dirname $(readlink -e $0))

chr=$1
ref=$2
contigs=$3
out_dir=$4

threads=8
if [ "$#" -gt 4 ]; then
    threads=$5
fi

mkdir -p $out_dir

align_prefix=$out_dir/align

minimap2 -a -H -t $threads -x asm20 -I 8G $ref $contigs | samtools view -b -F 4 > $align_prefix.bam

$root/../alignment_filter.py --min-len 4000 --min-idy 0.98 --use-all $align_prefix.bam 2> $align_prefix.cerr > $align_prefix.cout

$root/contig_info.py $contigs > $out_dir/contig.info

$root/gap_patch.py $out_dir/contig.info $align_prefix.cout 2> $out_dir/patch.log | sort | uniq | sed 's/^.*: //g' > $out_dir/patch.bed

#grep -Po "^.*caffold.*?\s" $out_dir/patch.bed | sort | uniq > $out_dir/relevant_ref.txt
cut -f 1 $out_dir/patch.bed | sort | uniq > $out_dir/relevant_ref.txt
seqtk subseq $ref $out_dir/relevant_ref.txt > $out_dir/for_patch.fasta
cat $contigs >> $out_dir/for_patch.fasta

echo ">$chr" > $out_dir/noformat.fasta
#bedtools getfasta -s -fi $out_dir/for_patch.fasta -bed $out_dir/patch.bed | grep -v ">" >> $out_dir/noformat.fasta
#TODO support streaming and get rid of intermediate fasta files
$root/merge_reads.py $out_dir/patch.bed $out_dir/for_patch.fasta | grep -v ">" | tr -d '\n' >> $out_dir/noformat.fasta
fold -c $out_dir/noformat.fasta > $out_dir/$chr.fasta
rm $out_dir/noformat.fasta

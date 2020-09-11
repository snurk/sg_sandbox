#!/bin/bash
set -e

if [ "$#" -lt 1 ]; then
    echo "script.sh <chromosome name>"
    echo "patch.bed file and cns.renamed.fasta files should be in the folder"
    exit 239
fi

chr=$1

grep -Po "^.*caffold.*?\s" patch.bed | sort | uniq > v07.txt
grep -Po "chrX_v0.7_0\s" patch.bed | sort | uniq >> v07.txt
seqtk subseq /data/PoreTenders/10x_v4CentromereXOnly/freebayes_serge/round2/splitBasedOnNCBI/new_contigs.fasta v07.txt > for_patch.fasta
cat cns.renamed.fasta >> for_patch.fasta

echo ">$chr" > noformat.fasta
#bedtools getfasta -s -fi for_patch.fasta -bed patch.bed | grep -v ">" >> noformat.fasta
~/git/ngs_scripts/gfakluge/contig_processing/merge_reads.py patch.bed for_patch.fasta | grep -v ">" >> noformat.fasta

~/git/ngs_scripts/contig_processing/contig_length_filter.py 1 noformat.fasta $chr.fasta

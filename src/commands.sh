#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "commands.sh <folder> <min_overlap> [weak_overlaps...]+"
    exit 239
fi

echo "Usage of file deprecated"
exit 2

#params: <folder> <min_ovl> <min_wea>
module load minimap2
module load samtools
module load quast

mkdir -p $1
cd $1

##~/git/ngs_scripts/canu_gap_analysis/homopolymer_compress.py ../defensin_region.fasta ../defensin_region.compressed.fasta
##~/git/ngs_scripts/canu_gap_analysis/homopolymer_compress.py ../canu.fasta canu.compressed.fasta

if [ -f ../reads.compressed.fasta.gz ]; then
    echo "File reads.compressed.fasta.gz was found in the folder and will be reused"
else
    ~/git/canu/Linux-amd64/bin/sqStoreDumpFASTQ -S ../assembly/asm.seqStore -nolibname -noreadname -fasta -compressed -trimmed -o ../reads.compressed.gz
fi

#rm -rf *.txt *.log *.gfa quast* ovl.paf *.sam *.nodes.fasta

~/git/ngs_scripts/gfakluge/pipe.sh ../reads.compressed.fasta.gz ../assembly $2 0.0001 ${@:3} &> pipe.log

#~/git/ngs_scripts/gfakluge/analyze_subgraph.sh . ../defensin_region.compressed.fasta

#~/git/ngs_scripts/gfakluge/node_coverage.sh simplified.subgraph.gfa ../reads.compressed.fasta.gz

#for g in simplified.gfa no_*.gfa processed.gfa ; do
#    name=$(basename $g .gfa)
#    awk '/^S/{print ">"$2"\n"$3}' $g | fold > $name.nodes.fasta
#    minimap2 -x asm10 -a -t 12 -o $name.nodes.sam ../defensin_region.compressed.fasta $name.nodes.fasta
#    ~/git/ngs_scripts/alignment_filter.py --min-len 50000 $name.nodes.sam 2> $name.filter.log | awk '{print $1}' > $name.nodes.tmp.txt
#    ~/git/ngs_scripts/alignment_filter.py --min-len 15000 --query-frac .90 $name.nodes.sam 2>> $name.filter.log | awk '{print $1}' >> $name.nodes.tmp.txt
#    sort $name.nodes.tmp.txt | uniq > $name.nodes.txt
#    rm $name.nodes.tmp.txt
#    ~/git/ngs_scripts/gfakluge/neighborhood $g $name.subgraph.gfa $name.nodes.txt 5 &> $name.nodes_rad_5.log
#done
#
##awk '/^S/{print ">"$2"\n"$3}' simplified.gfa | fold > simplified.nodes.fasta
##
##minimap2 -x asm10 -a -t 12 -o simplified.nodes.sam ../defensin_region.compressed.fasta simplified.nodes.fasta
##
###FIXME tweak the filter. At least 0.9
##~/git/ngs_scripts/alignment_filter.py simplified.nodes.sam 5000 .80 simplified.nodes.filtered.sam 1000 2> simplified.filter.log | awk '{print $1}' | sort | uniq > simplified.nodes.txt
##
##~/git/ngs_scripts/gfakluge/neighborhood simplified.gfa simplified.subgraph.gfa simplified.nodes.txt 5 &> simplified.nodes_rad_5.log
#
##quast.py -o quast_cmp -R ../defensin_region.compressed.fasta -t 12 --min-identity 99 --min-alignment 1000 --ambiguity-usage all simplified.nodes.fasta ../canu.compressed.fasta #../flye.utg.compressed.fasta ../flye.compressed.fasta
#
##minimap2 -x asm10 -a -I 10G -t 12 simplified.nodes.fasta ../reads.compressed.fasta.gz | samtools view -q 20 > alignment.sam
##
##~/git/ngs_scripts/alignment_filter.py --min-len 5000 --query-frac 0.9 --filtered alignment.filtered.sam alignment.sam &> read.filter.log
##
##~/git/ngs_scripts/alignment_stats/compute_coverage.py alignment.filtered.sam &> alignment.filtered.cov.info
#
#awk '/^S/{print ">"$2"\n"$3}' simplified.subgraph.gfa | fold > simplified.subgraph.nodes.fasta
#
#minimap2 -x asm10 -a -I 10G -t 12 simplified.subgraph.nodes.fasta ../reads.compressed.fasta.gz | samtools view -q 20 -b > subgraph.alignment.bam
#
#~/git/ngs_scripts/alignment_filter.py --min-len 5000 --query-frac 0.9 --filtered subgraph.alignment.filtered.bam subgraph.alignment.bam &> read.filter.log
#
#~/git/ngs_scripts/alignment_stats/compute_coverage.py subgraph.alignment.filtered.bam &> subgraph.alignment.filtered.cov.info

cd -

#!/bin/bash
set -e

####Launch consensus#########
mkdir -p consensus
cd consensus

#Uncomment for consensus over all unitigs
#ln -s ../resolved_mapping.txt backbone_layout.txt

sbatch --ntasks 1 --mem 20G --cpus-per-task 10 --time 6:00:00  ~/git/ngs_scripts/gfakluge/consensus/consensus.sh

####Combine contigs#########
conda activate

cat cns-*.cns.fasta > cns.result.fasta

awk '{print $1,NR}' backbone_layout.txt > cns_name_map.txt

~/git/ngs_scripts/canu_gap_analysis/rename_reads.py cns.result.fasta cns_name_map.txt cns.renamed.fasta

~/git/ngs_scripts/contig_processing/contig_info.py cns.renamed.fasta cns.renamed.info

####Alignments for gap closing##########
CHR=chr4

v07_root=/data/PoreTenders/10x_v4CentromereXOnly/freebayes_serge/round2/splitBasedOnNCBI

cat $v07_root/chm13_v06_chr_contig_assignments.txt | awk -v C=$CHR '{V="chr"$2; if (V == C) print $1}' > v07.relevant.txt

seqtk subseq $v07_root/new_contigs.fasta v07.relevant.txt > v07.relevant.fasta

#module load minimap2
#minimap2 -ax asm20 -H -t 16 cns.renamed.fasta > v07.sam

sbatch --ntasks 1 --mem 20G --cpus-per-task 6 --time 4:00:00 ~/git/ngs_scripts/gfakluge/contig_processing/align.sh cns.renamed.fasta v07.relevant.fasta v07.bam 5

#sbatch --ntasks 1 --mem 20G --cpus-per-task 6 --time 4:00:00 ~/git/ngs_scripts/gfakluge/contig_processing/align.sh cns.renamed.fasta $v07_root/new_contigs.fasta v07.full.bam 5

~/git/ngs_scripts/alignment_filter.py --min-len 4000 --min-idy 0.98 --use-all v07.bam 2> v07_align.cerr > v07_align.cout

~/git/ngs_scripts/gfakluge/contig_processing/gap_patch.sh $CHR cns.renamed.info v07_align.cout test_patch.bed

#/data/korens/devel/utils/bamparse/samToErrorRate v07.sam $CHR.v07.fasta 200 | awk '{if (match($1, "#")) { print $0; } else { if ($9 == 0) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t1\t"$12-$11"\t"$12-$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19}}' | sort -nk10 > v07.paf
#
#awk '{print "Alignment of query",$1,"( length",$8,") : [",$6,"-",$7,"], orientation:",$9,"to",$2,"( length",$12,") : [",$10,"-",$11,"]. Identity --",$4}' v07.paf > v07.paf.readable

#grep "^chr19.3" $CHR.paf | awk '{print "Alignment of query",$1,"( length",$8,") : [",$6,"-",$7,"], orientation:",$9,"to",$2,"( length",$12,") : [",$10,"-",$11,"]. Identity --",$4}' | grep Super-Scaffold_466_0

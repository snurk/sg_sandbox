#!/bin/bash

chr=$1
out=$2
module load samtools
module load minimap2

if [ ! -s $out/merged.sorted.bam.bai ]; then
   samtools merge -O bam $out/merged.sorted.bam $out/map/*.bam
   samtools index $out/merged.sorted.bam
fi

if [ ! -s $out/reads.cut ]; then
   samtools view -h $out/merged.sorted.bam | k8 `which paftools.js` sam2paf -p - |awk '{if ($10/$11 >= 0.85 && $4-$3 >= 10000) print $0}' | grep -w $chr > $out/merged.sorted.paf

   cat $out/merged.sorted.paf |awk '{print $1"\t0\t"$2}'|sort |uniq > $out/reads.cut
   EXP=`wc -l $out/reads.cut |awk '{print $1}'`

   > $out/reads.fastq
   for i in `cat input.fofn`; do
      seqtk subseq $i $out/reads.cut 2> $out/reads_extract.err >> $out/reads.fastq
   done
   pigz $out/reads.fastq
   ACT=`seqtk comp $out/reads.fastq.gz |wc -l`

   if [ $ACT -ne $EXP ]; then
      echo "Error: expected $EXP reads but only extracted $ACT"
      rm $out/reads.*
      exit
   fi
   zcat $out/reads.fastq.gz |awk '{if ((NR-1)%4 <= 1) print $0}' |sed -r s/'A{1,}+'/A/g |sed -r s/'C{1,}+'/C/g |sed -r s/'G{1,}+'/G/g |sed -r s/'T{1,}+'/T/g |tr '@' '>' | awk '{print $1}' | fold -c | pigz -c - > $out/reads_compressed.fasta.gz
   ACT=`seqtk comp $out/reads_compressed.fasta.gz| wc -l`
   if [ $ACT -ne $EXP ]; then
      echo "Error: expected $EXP reads but only extracted $ACT"
      rm $out/reads.*
      exit
   fi
fi

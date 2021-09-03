#!/bin/bash

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <layouts> <graph in gfa> <file with fasta of nodes> <file with graphaligner alignments to gfa> <ont reads fastq>"
    exit 1
fi

base_path=$(dirname $(readlink -e $0))

layouts=$1
gfa=$2
fasta=$3
gaf=$4
reads=$5

# get the sequence for gaps and make candidate read list
echo "Getting sequences"
$base_path/../join_segments.py $fasta $layouts $gfa allpaths.fasta
if [ ! -e allpaths.fasta ]; then
   echo "Error: could not generate candidate fasta, check for errors above and confirm your traversals"
   exit
fi
cat allpaths.fasta | awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta"); } print $0 > filename }'
rm -f allpaths.fasta

echo "Getting read lists"
for gap in `cat $layouts |awk '{print substr($1, 1, index($1, "v")-1)}'|sort |uniq`; do
echo "Processing gap $gap"
   cat $layouts |grep $gap | awk '{print $NF}' |sed s/+//g |sed s/-//g |tr ',' '\n' |sort |uniq > $gap.nodes
   grep -w -f $gap.nodes $gaf | awk '{gsub(">", ",", $2); gsub("<", ",", $2); if (gsub(",", "", $2) > 1) { print $1"\t0\t10000000"}}' |sort |uniq > $gap.reads.cut
done

for i in `ls *.cut`; do
   prefix=`echo $i | sed s/.reads.cut//g |sed s/.cut//g`
   if [ ! -e $prefix.fastq.gz ]; then
echo "Need to extract reads for gap $prefix"
      seqtk subseq $reads $i 2> $prefix.err |pigz -c - > $prefix.fastq.gz
   else
      echo "Already extracted"
      EXP=`wc -l $i |awk '{print $1}'`
      NUM=`seqtk comp $prefix.fastq.gz |wc -l |awk '{print $1}'`
      if [ $EXP -ne $NUM ]; then
         echo "Error: extracted set for $i has only $NUM reads, expected $EXP"
         exit
      fi
   fi
done

# compress the reads
for i in `ls *.fastq.gz`; do
   prefix=`echo $i |sed s/.fastq.gz//g`
   if [ ! -e ${prefix}_compress.fasta.gz ];then
      echo "Compressing $i"
      zcat $i  |awk '{if ((NR-1)%4 <= 1) print $0}' |sed -r s/'A{1,}+'/A/g |sed -r s/'C{1,}+'/C/g |sed -r s/'G{1,}+'/G/g |sed -r s/'T{1,}+'/T/g |tr '@' '>' |pigz -c - >${prefix}_compress.fasta.gz
   fi
done

for i in `ls *gap*fasta|grep -v compress`; do
   name=`echo $i |sed s/.fasta//g`
   prefix=`echo $name |awk -F "v" '{print $1}'`
   if [ ! -e $name.paf ]; then
echo "Mapping $name $prefix.gz"
   /data/korens/devel/Winnowmap/bin/winnowmap -z150 -t 16 -ax map-ont -H -k 19 $name.fasta ${prefix}_compress.fasta.gz > $name.sam
   /data/korens/devel/utils/bamparse/samToErrorRate $name.sam $name.fasta 1000 | awk '{if (match($1, "#")) { print $0; } else { if ($9 == 0) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t1\t"$12-$11"\t"$12-$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19}}' |sort -nk10|awk '{if ($4 > 50) print $0}' > $name.paf
   fi
done

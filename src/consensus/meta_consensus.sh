#!/bin/bash
set -e

bin=~/git/canu2/Linux-amd64/bin

#seqstore=../../../assembly/asm.seqStore

seqstore=../../assembly/asm.seqStore
ovlstore=../../assembly/unitigging/asm.ovlStore

#if [ ! -f backbone_layout.txt ] ; then
#  ~/git/ngs_scripts/gfakluge/consensus/resolve_layouts.sh layout.txt ../../8k_upd/resolved_mapping.txt
#fi

cnspre=./cns

$bin/layoutReads -S $seqstore -O $ovlstore -eg 0.001 -eM 0.001 -R ../resolved_mapping.txt -gs 300000000 -o $cnspre &> layoutReads.log

#  Magic values
#    partitionSize    - make partitions this big, relative to the biggest contig
#                       0.10 - each big contig gets its own partition
#    partitionScaling - expect contigs to expand by this much due to homopoly compression
#                       1.50 - 50% expansion of contig length, affects memory estimate
#    partitionTigs    - put this many reads (as fraction of total) in a partition
#                       0.01 - puts small contigs into many partitions

errorRate=0.05

if [ ! -e ${cnspre}.ctgStore/partitioning ] ; then
  $bin/utgcns -V \
    -S ${seqstore} \
    -T ${cnspre}.ctgStore 1 \
    -partition 0.10 1.50 0.01 \
  &> ${cnspre}.ctgStore/partitioning.log
fi

parts=`cd ${cnspre}.ctgStore ; ls partition.???? | sed s/partition.//`
for pp in $parts ; do
  echo "#!/bin/bash" > ./${cnspre}_launch.${pp}.sh
  echo \
  $bin/utgcns \
    -R ${cnspre}.ctgStore/partition.${pp} \
    -T ${cnspre}.ctgStore 1 \
    -P ${pp} \
    -O ${cnspre}-${pp}.cns \
    -A ${cnspre}-${pp}.cns.fasta \
    -maxcoverage 50 \
    -e ${errorRate} \
    -threads 8 >> ./${cnspre}_launch.${pp}.sh

  sbatch --ntasks 1 --mem 100G --cpus-per-task 10 --time 12:00:00 ./${cnspre}_launch.${pp}.sh
done



#!/bin/bash
set -eou

if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <canu bin> <backbone layout> <seqstore> <ovlstore> <out folder>"
    exit 1
fi

bin=$(readlink -e $1)
backbone_layout=$(readlink -e $2)
seqstore=$(readlink -e $3)
ovlstore=$(readlink -e $4)

mkdir -p $5
cd $5

#bin=~/git/canu2/build/bin
#
##seqstore=../../../assembly/asm.seqStore
#
#seqstore=../../raw_seqstore_works/asm.seqStore
##ovlstore=../../assembly/unitigging/asm.ovlStore
#ovlstore=../../assembly/unitigging/4-unitigger/asm.0.all.ovlStore
#
#if [ ! -f backbone_layout.txt ] ; then
#  $(dirname $0)/resolve_layouts.sh layout.txt ../resolved_mapping.txt
#fi

cnspre=./cns

$bin/layoutReads -S $seqstore -O $ovlstore -eg 0.00001 -eM 0.00001 -R $backbone_layout -seed 239 -gs 200000000 -o $cnspre &> layoutReads.log

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
  $bin/utgcns -V -V -V \
    -R ${cnspre}.ctgStore/partition.${pp} \
    -T ${cnspre}.ctgStore 1 \
    -P ${pp} \
    -O ${cnspre}-${pp}.cns \
    -A ${cnspre}-${pp}.cns.wip.fasta \
    -maxcoverage 50 \
    -e ${errorRate} \
    -threads 8 "&>" ${cnspre}-${pp}.cns.log >> ./${cnspre}_launch.${pp}.sh

  echo "mv ${cnspre}-${pp}.cns.wip.fasta ${cnspre}-${pp}.cns.fasta" >> ./${cnspre}_launch.${pp}.sh

  sbatch --ntasks 1 --mem 100G --cpus-per-task 10 --time 24:00:00 ./${cnspre}_launch.${pp}.sh
done

cd -

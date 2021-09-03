#!/bin/bash
EXE="/data/Phillippy/tools/Winnowmap/bin/winnowmap"
OPT="map-ont"
REF=$1
JOBID=$SLURM_ARRAY_TASK_ID
THREADS=`echo "$SLURM_CPUS_PER_TASK" |awk '{print int($1/2)}'`

QRY=`head -n $JOBID input.fofn |tail -n 1`
PREFIX=`basename $QRY |sed s/.fastq.gz//g |sed s/.fq.gz//g`
PREFIX="$2/$PREFIX"

module load samtools

if [ ! -e ${PREFIX}.bam ]; then
   /usr/bin/time ${EXE} --MD -H -ax $OPT -t $THREADS $REF $QRY > ${PREFIX}.sam
   /usr/bin/time samtools sort -m 20G -@4 -O BAM ${PREFIX}.sam > ${PREFIX}.bam
fi

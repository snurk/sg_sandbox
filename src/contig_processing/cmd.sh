#!/bin/bash
set -e

# PREPARE

hicanu=/data/Phillippy/projects/chm13/T2T_Finishing_Workshop/canu_master0424/asm.contigs.fasta

mkdir -p alignments

for f in ~/globus/assemblies/drafts/20200611/*.fasta.gz ~/globus/assemblies/drafts/20200602/*.fasta.gz /gpfs/gsfs7/users/nurks2/gfa_works/chm13/new_oea/chr_reconstruction_8k/*/*v1p5.fasta.gz ~/data/gfa_works/chm13/symm_oea/chr_reconstruction_8k/chr5/*.fasta.gz /gpfs/gsfs7/users/nurks2/gfa_works/chm13/symm_oea/chr_reconstruction_8k/telomere_fix/chr6_centroFlye/chr6_centroFlye.fasta.gz ; do
    chr=$(basename $f .fasta.gz)
    if [ "$chr" = "chr6" ]; then
        continue
    fi
    if [ "$chr" = "chr17" ]; then
        continue
    fi
    if [ "$chr" = "chr8" ]; then
        continue
    fi
    if [ "$chr" = "chrX" ]; then
        continue
    fi
    if [ "$chr" = "chr8.glennis" ]; then
        continue
    fi

    echo "Processing chromosome $chr from $f"

    if [ ! -f alignments/${chr}_to_hicanu.bam ] ; then
        sbatch --ntasks 1 --mem 80G --cpus-per-task 5 --time 12:00:00 ./align.sh $f $hicanu alignments/${chr}_to_hicanu.bam 5
    fi

    mkdir -p $chr
    if [ ! -f $chr/init.info ] ; then
        ~/git/ngs_scripts/contig_processing/contig_info.py $f > $chr/init.info
    fi

    rm -f $chr/unpatched.fasta.gz
    ln -s $(readlink -e $f) $chr/unpatched.fasta.gz

    #if [ ! -f alignments/hicanu_to_$chr.bam ] ; then
    #    sbatch --ntasks 1 --mem 80G --cpus-per-task 20 --time 12:00:00 ./align.sh $hicanu $f alignments/hicanu_to_$chr.bam 20
    #fi
done

# PATCH

hicanu=/data/Phillippy/projects/chm13/T2T_Finishing_Workshop/canu_master0424/asm.contigs.fasta

rm -f for_patch.fasta*
cat $hicanu > for_patch.fasta
zcat */unpatched.fasta.gz >> for_patch.fasta

for f in alignments/*_to_hicanu.bam ; do
    chr=$(basename $f _to_hicanu.bam)
    #if [ "$chr" = "chr6" ]; then
    #    continue
    #fi
    #if [ "$chr" = "chr17" ]; then
    #    continue
    #fi
    #if [ "$chr" = "chr8.glennis" ]; then
    #    continue
    #fi
    #if [ "$chr" = "chr8" ]; then
    #    continue
    #fi
    echo "Processing chromosome $chr"
    #Diagnostic output
    #samtools view -h alignments/${chr}_to_hicanu.bam > alignments/${chr}_to_hicanu.sam
    #/data/korens/devel/utils/bamparse/samToErrorRate alignments/${chr}_to_hicanu.sam $hicanu 200 | grep -v "#Large" | awk '{if (match($1, "#")) { print $0; } else { if ($9 == 0) print $0; else print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t1\t"$12-$11"\t"$12-$10"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19}}' | sort -nk10 > ${chr}/to_hicanu.paf
    #awk '{print "Alignment of query",$1,"( length",$8,") : [",$6,"-",$7,"], orientation:",$9,"to",$2,"( length",$12,") : [",$10,"-",$11,"]. Identity --",$4}' ${chr}/to_hicanu.paf > ${chr}/to_hicanu.readable

    ~/git/ngs_scripts/alignment_filter.py --min-idy 0.99 --min-len 10000 --use-all alignments/${chr}_to_hicanu.bam 2> ${chr}/parse.log > ${chr}/to_hicanu.align

    ./telomere_patch.py ${chr}/to_hicanu.align > ${chr}/telomere_patch.bed
done &> patch_gen.log

sed -i 's/tig00000001/chr11/g' chr11/telomere_patch.bed

conda activate
module load bedtools
for f in */telomere_patch.bed ; do
   name=$(basename $(dirname $f) v1p5)
   name=$(basename $name _centroFlye)
   if [ -s $f ]; then
       echo "chr $name WILL be patched";
       cd $(dirname $f)
       ../patch_from_bed.sh ${name}v2
       cd -
   else
       echo "chr $name WILL NOT be patched"
   fi
done

mkdir -p patches
cd patches
for f in ../*/telomere_patch.bed ; do
   name=$(basename $(dirname $f) v1p5)
   name=$(basename $name _centroFlye)
   if [ -s $f ]; then
       echo "chr $name WILL be patched";
       ln -s $(dirname $f)/telomere_patch.bed ${name}.telomere.bed
   else
       echo "chr $name WILL NOT be patched"
   fi
done
cd -

gzip */*v2.fasta

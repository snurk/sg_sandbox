#!/bin/bash
set -e

module load quast

g=simplified.subgraph.gfa
name=$(basename $g .gfa)

cd fragmented_region

awk '/^S/{print ">"$2"\n"$3}' $g | fold > $name.nodes.fasta
quast.py -t 12 --large --min-identity 98. --min-alignment 5000 -R /home/nurks2/data/fragmentation_invest/chm13/analysis/fragmented_region_chrX/compressed.region.fasta -o ${name}.quast $name.nodes.fasta /home/nurks2/data/fragmentation_invest/chm13/analysis/fragmented_region_chrX/tigs.compressed.fasta

cd -

cd defensin

awk '/^S/{print ">"$2"\n"$3}' $g | fold > $name.nodes.fasta
quast.py -t 12 --large --min-identity 98. --min-alignment 5000 -R ~/data/gfa_works/chm13/extended_region.v3.compressed.fasta -o ${name}.quast $name.nodes.fasta ~/data/gfa_works/chm13/hicanu_defensin.compressed.fasta

cd -

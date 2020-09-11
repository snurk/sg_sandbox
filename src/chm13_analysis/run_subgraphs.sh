#!/bin/bash
set -e

module load minimap2
module load samtools

mkdir -p fragmented_region
cd fragmented_region
~/git/ngs_scripts/gfakluge/analyze_subgraph.sh ../simplified.gfa ~/data/fragmentation_invest/chm13/analysis/fragmented_region_chrX/compressed.region.fasta
cd -

mkdir -p defensin
cd defensin
~/git/ngs_scripts/gfakluge/analyze_subgraph.sh ../simplified.gfa ~/data/gfa_works/chm13/extended_region.v3.compressed.fasta
cd -

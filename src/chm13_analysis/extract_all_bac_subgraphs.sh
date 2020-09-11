#!/bin/bash
set -e

module load minimap2
module load samtools

mkdir -p all_bacs
cd all_bacs
~/git/ngs_scripts/gfakluge/analyze_bac_subgraph.sh ../simplified.gfa ../../bacs.compressed.fasta
cd -

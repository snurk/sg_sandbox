#!/bin/bash
set -e

cat ../../../8k/chr_analysis/chrY/final_compress_mapping.txt ../../../8k/mapping.txt > old_mapping.txt

cat to_lift.txt old_mapping.txt > tmp_mapping.txt
awk '{print $1}' to_lift.txt > tmp_names.txt
~/git/ngs_scripts/gfakluge/consensus/resolve_layouts.py --miniasm ../../../8k/microasm.gfa tmp_names.txt tmp_mapping.txt > old_backbone_layout.txt
~/git/ngs_scripts/gfakluge/coord_uplift.py old_backbone_layout.txt ../../resolved_mapping.txt #> uplift.log

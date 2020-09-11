#!/bin/bash
set -e

for g in simplified*000.gfa simplified.final_*.gfa ; do
    name=$(basename $g .gfa)
    d=$name

    rm -rf $d
    mkdir -p $d
    ~/git/ngs_scripts/bogart_gfa/resolve_mapping.py $g microasm.gfa mapping.txt > $d/resolved_mapping.txt

    ~/git/ngs_scripts/gfakluge/assign_coverage.py $d/resolved_mapping.txt min_read.cov > $d/simplified.cov

    #awk '/^S/{print ">"$2"\n"$3}' simplified.gfa | fold > simplified.nodes.fasta

    cp $g $d/no_cov.gfa
    ~/git/ngs_scripts/gfakluge/inject_coverage.py $d/no_cov.gfa $d/simplified.cov > $d/simplified.gfa
    rm $d/no_cov.gfa

    echo -e "H\tVN:Z:1.0" > $d/simplified.noseq.gfa
    ~/git/gfatools/gfatools view -S $d/simplified.gfa >> $d/simplified.noseq.gfa
done

#!/bin/bash
set -e

scripts_root=$(dirname $(readlink -e $0))/..

for g in simplified*000.gfa simplified.final_*.gfa ; do
    name=$(basename $g .gfa)
    d=$name

    rm -rf $d
    mkdir -p $d
    $scripts_root/resolve_layouts.py $g mapping.txt --miniasm microasm.gfa > $d/resolved_mapping.txt

    $scripts_root/assign_coverage.py $d/resolved_mapping.txt min_read.cov > $d/simplified.cov

    #awk '/^S/{print ">"$2"\n"$3}' simplified.gfa | fold > simplified.nodes.fasta

    cp $g $d/no_cov.gfa
    $scripts_root/inject_coverage.py $d/no_cov.gfa $d/simplified.cov > $d/simplified.gfa
    rm $d/no_cov.gfa

    echo -e "H\tVN:Z:1.0" > $d/simplified.noseq.gfa
    $scripts_root/../gfacpp/gfatools/gfatools view -S $d/simplified.gfa >> $d/simplified.noseq.gfa
done

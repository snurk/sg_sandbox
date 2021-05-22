#!/bin/bash
set -e

#scripts_root=$(dirname $(readlink -e $0))/..
scripts_root=~/git/sg_sandbox/src

#for g in simplified*000.gfa simplified.final_*.gfa ; do
for g in simplified*000.gfa ; do
    name=$(basename $g .gfa)
    d=$name

    rm -rf $d
    mkdir -p $d
    $scripts_root/resolve_layouts.py $g mapping.txt --miniasm microasm.gfa > $d/resolved_mapping.txt

    $scripts_root/assign_coverage.py $d/resolved_mapping.txt min_read.cov > $d/simplified.cov

    $scripts_root/inject_coverage.py $g $d/simplified.cov > $d/simplified.gfa

    awk '/^S/{print ">"$2"\n"$3}' $d/simplified.gfa | fold > $d/simplified.nodes.fasta

    echo -e "H\tVN:Z:1.0" > $d/simplified.noseq.gfa
    $scripts_root/../gfacpp/gfatools/gfatools view -S $d/simplified.gfa >> $d/simplified.noseq.gfa

    #cd $d
    #$scripts_root/chm13_analysis/mashmap.sh simplified.nodes.fasta &> mashmap.log
    #$scripts_root/chm13_analysis/chr_analysis.sh &> chr_analysis.log
    #cd -
done

#!/bin/bash
set -eou

#if [ "$#" -lt 4 ]; then
#    echo "Usage: $(basename $0) <processed.gfa> <init.gfa> <node_composition.txt> <read_coverage.txt>"
#    exit 1
#fi

awk '/^S/{print ">"$2"\n"$3}' simplified.gfa | fold > simplified.nodes.fasta

echo -e "H\tVN:Z:1.0" > simplified.noseq.gfa
~/git/sg_sandbox/gfacpp/gfatools/gfatools view -S simplified.gfa >> simplified.noseq.gfa

~/git/sg_sandbox/src/chm13_analysis/mashmap.sh simplified.nodes.fasta &> mashmap.log
~/git/sg_sandbox/src/chm13_analysis/chr_analysis.sh &> chr_analysis.log

#!/bin/bash
set -eou

root=$(dirname $(readlink -e $0))/..

n=$1
cd $n
$root/join_segments.py ../simplified.nodes.fasta layout.txt ../simplified.noseq.gfa joined.fasta &> join.log

cat ../$n\.*/cns.renamed.fasta > stitch.fasta

$root/contig_processing/contig_info.py stitch.fasta stitch.info
$root/contig_processing/contig_info.py joined.fasta joined.info

$root/consensus/join_ctgs_parasail.py stitch.fasta stitch_plan.txt stitched.fasta 30000 300 &> stitch.log
$root/contig_processing/contig_info.py stitched.fasta stitched.info

ln -s stitched.fasta cns.renamed.fasta
cd -

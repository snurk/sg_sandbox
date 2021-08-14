#!/bin/bash
set -e

if [ "$#" -lt 1 ]; then
    echo "script.sh <flank size> <penalty>"
    exit 239
fi

~/git/ngs_scripts/join_ctgs.py stitch.fasta stitch_plan.txt stitched.fasta $1 $2 &> stitch.log

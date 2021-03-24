#!/bin/bash

module load mashmap

mashmap -s 100000 -t 16 --pi 95 -q $1 -r ~/data/gfa_works/hg38.chronly.compressed.fasta -o $2

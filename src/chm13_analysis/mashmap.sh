#!/bin/bash
set -eou

module load mashmap

mashmap -f map -s 10000 -t 16 --pi 95 -q $1 -r $2 -o $3

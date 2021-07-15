#!/bin/bash

module load mashmap

mashmap -s 100000 -t 16 --pi 95 -q $1 -r $2 -o $3

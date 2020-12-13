#!/bin/bash
set -eou

grep "^S" | awk '{print $2,$5}' | sed 's/ll:f://g' | sed 's/\.*//g'

#!/bin/bash
#set -eou

grep "^S" | sed 's/^S\s//g' | sed 's/\s.*ll:f:/ /g' | awk '{print $1,$2}'

#!/bin/bash
#set -eou

grep "^S" | awk '{print $2,$4,$5}' | sed 's/..:i://g' | awk '{print $1,$3/$2}'

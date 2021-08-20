#!/bin/bash
set -e

orig=$1
oea=$2

bin=/home/nurks2/git/canu2/build/bin/

$bin/sqStoreDumpFASTQ -S $orig/asm.seqStore -nolibname -fasta -normal -o reads.gz

cp $oea/asm.seqStore/readNames.txt oeaNames.txt

awk '{print $2}' oeaNames.txt > names.txt

seqtk subseq reads.fasta.gz names.txt | gzip --fast > filtered.fasta.gz

$bin/sqStoreCreate \
  -o asm.seqStore \
  -minlength 1000 \
  -genomesize 3100000000 \
  -coverage   100000 \
  -bias       0 \
  -raw -trimmed -pacbio-hifi all_data filtered.fasta.gz

touch asm.seqStore/homopolymerCompression

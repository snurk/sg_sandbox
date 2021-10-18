#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
  echo "Usage: $0 <raw reads FASTA (.gz)> <second round seqstore> <out prefix>"
    exit 1
fi

raw_reads=$1
second_seqstore=$2
prefix=$3

mkdir -p $(dirname $prefix)

#$CANU_BIN/sqStoreDumpFASTQ -S $first_seqstore -nolibname -fasta -normal -o $prefix.reads.gz

awk '{print $2}' $second_seqstore/readNames.txt > $prefix.names.txt

seqtk subseq $raw_reads $prefix.names.txt | gzip --fast > $prefix.filtered.fasta.gz

rm -rf $prefix
rm -rf $prefix.BUILDING

$CANU_BIN/sqStoreCreate \
  -o $prefix.BUILDING \
  -minlength 0 \
  -genomesize 3100000000 \
  -coverage   100000 \
  -bias       0 \
  -raw -trimmed -pacbio-hifi all_data $prefix.filtered.fasta.gz \
&& mv $prefix.BUILDING $prefix

echo "Looking at difference between $second_seqstore/readNames.txt and $prefix/readNames.txt"
diff -q <(sed 's/ .*//g' $second_seqstore/readNames.txt) <(sed 's/ .*//g' $prefix/readNames.txt)

rm -f $prefix.filtered.fasta.gz $prefix.names.txt

touch $prefix/homopolymerCompression

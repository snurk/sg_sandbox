#!/usr/bin/env python
import sys
import gzip
from Bio import SeqIO

# the patch file is expected to be 0-based and not including on the last one ala BED

def reverse_complement(kmer):
   complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a' : 't', 'c' : 'g', 'g' : 'c', 't' : 'a'}
   return str("".join(complement.get(base, base) for base in reversed(kmer)))

if len(sys.argv) < 3:
    print("Usage: %s <bed file> <fasta>" % sys.argv[0], file=sys.stderr)
    sys.exit(1)

# the first file tells us all the order
ids = []

with open('%s'% (sys.argv[1])) as f:
   for line in f:
      if line.startswith("#"):
         continue

      l = line.strip().split()
      if len(l) > 3 and l[5] == '-':
         ids.append("%s %s %s"%(l[0], int(l[2]), int(l[1])))
      else:
         ids.append("%s %s %s"%(l[0], int(l[1]), int(l[2])))

# this file gives us sequences, load into ram
seqs = {}
recs = [(rec.name, str(rec.seq)) for rec in SeqIO.parse(open(sys.argv[2]), "fasta")]
for name, seq in recs:
   seqs[name] = seq
seqs["gap"] = "N" * 10000000

print(">merged")
for i in ids:
   i = i.split()
   start = int(i[1])
   end = int(i[2])
   seq = seqs[i[0]]

   if start < end:
      seq = seqs[i[0]][start:end]
   else:
      seq = reverse_complement(seqs[i[0]][end:start])
   print(seq)

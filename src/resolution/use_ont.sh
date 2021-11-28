#!/bin/bash
set -eou

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <ONT reads> <graph in gfa> <file with names of single-copy unitigs> <out folder> [alignment threads, default = 8]"
    exit 1
fi

threads=8
if [ "$#" -gt 4 ]; then
    threads=$5
fi

base_path=$(dirname $(readlink -e $0))

reads=$1
g=$2
unique=$3
out=$4

mkdir -p $out

if [ ! -f $out/ga_ont.gaf ] ; then
    echo "Aligning reads with GraphAligner"
    #sbatch --ntasks 1 -W --mem 80G --cpus-per-task 24 --time 24:00:00 $base_path/ga_ont.sh $g $reads $out/ga_ont $threads
    $base_path/ga_ont.sh $g $reads $out/ga_ont $threads
fi

echo "Proceeding with resolution"
awk '$1!="S"{print;}$1=="S"{print "S\t" $2 "\t*\tLN:i:" length($3) "\t" $4 "\t" $5 "\t" $6 "\t" $7;}' < $g > $out/noseq.gfa

#grep "^S" $out/noseq.gfa | sed 's/LN:i://g' | sed 's/RC:i://g' | awk '{print $2,$4,$5/$4}' > nodecovs.csv
#awk '{if (($2 >= 20000 && $3 >= 15 && $3 <= 30) || ($2 >= 50000 && $3 >= 10 && $3 <= 35) || $2 >= 100000) print $1}' < nodecovs.csv > unique_nodes.txt

sed 's/ <unknown description>//g' $out/ga_ont.gaf | awk -F "\t" '{split($1, name, " "); print name[1],($4-$3)/$2,$16,$6,$8,$9}' | sed 's/id:f://g' | awk -v OFS='\t' '{if ($2 > 0.9 && $3 >= 0.9) print $1,$4,$5,$6}' > $out/ont_filtered.tsv

$base_path/post_process_gaf.py --trusted-overhang 5000 --gaf-paths $out/ont_filtered.tsv $out/noseq.gfa $out/ont_processed.tsv &> $out/ont_process.log

#~/git/ngs_scripts/post_process_gaf.py --trusted-overhang 5000 ont_filtered.tsv simplified.noseq.gfa ont_processed.format.tsv

mikko_scripts=$base_path/../../tangle-resolution/scripts
cut -f 2 < $out/ont_processed.tsv | $mikko_scripts/find_bridges.py $unique > $out/ont_bridges.txt

grep -v '(' < $out/ont_bridges.txt | grep -vP '^$' | $mikko_scripts/connections_filter.py | sort > $out/bridging_seq_all.txt
#grep -v '(' < $out/ont_bridges.txt | grep -vP '^$' | $mikko_scripts/remove_wrong_connections_2.py | sort > $out/bridging_seq_all.txt


$mikko_scripts/pick_majority_bridge.py < $out/bridging_seq_all.txt > $out/bridging_seq_picked.txt
# $mikko_scripts/get_clean_connections.py bridging_seq_picked.txt | less

$mikko_scripts/forbid_unbridged_tangles.py $unique $out/noseq.gfa $out/bridging_seq_picked.txt > $out/forbidden_ends.txt
#$mikko_scripts/connect_uniques.py $out/noseq.gfa $out/forbidden_ends.txt $out/bridging_seq_picked.txt > $out/resolved.gfa 2> $out/connect_uniques.log
$mikko_scripts/connect_uniques.py $g $out/forbidden_ends.txt $out/bridging_seq_picked.txt > $out/resolved.gfa 2> $out/connect_uniques.log

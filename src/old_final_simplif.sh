#!/bin/bash
set -eou

if [ "$#" -lt 5 ]; then
    echo "script.sh <gfa> <utg mapping> <utg reads (.gfa with a-lines)> <read cov> <bubble_diff>"
    exit 239
fi

gfa=$1
utg_mapping=$2
utg_reads=$3
read_cov=$4
#bubble_diff=2000
bubble_diff=$5

scripts_root=$(dirname $(readlink -e $0))
algo_root=$scripts_root/../gfacpp/build

echo "Processing graph $gfa"
echo "Bubble diff length set to $bubble_diff"

grep -v "^a" $gfa | grep -v "^x" > simplified.wip.gfa

cp $utg_mapping mapping.txt

#Compaction rounds counter
cnt=100

final_it_cnt=3

for final_it in $(seq 1 $final_it_cnt) ; do

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Killing loops"
    $algo_root/loop_killer simplified.wip.gfa uncompressed.gfa simplified.wip.cov 30 > loop_killer_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Removing simple bulges ${final_it}"
    $algo_root/simple_bulge_removal simplified.wip.gfa uncompressed.gfa -l 10000 -d 3 --min-alt-ovl 10000 &> simple_br_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Removing shortcuts"
    $algo_root/shortcut_remover simplified.wip.gfa uncompressed.gfa simplified.wip.cov 30 8 &> shortcut_removal_${final_it}.log
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Final bubble removal ${final_it}"
    $algo_root/bubble_removal simplified.wip.gfa uncompressed.gfa 60000 $bubble_diff &> final_bfinder_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Removing low frequency hets in unique areas"
    $algo_root/simple_bulge_removal simplified.wip.gfa uncompressed.gfa --coverage simplified.wip.cov \
        --min-alt-ovl 10000 --max-unique-cov 30 --max-cov-ratio 0.33 -l 30000 -d 5000 &> low_freq_br_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Removing other variants in unique areas"
    $algo_root/simple_bulge_removal simplified.wip.gfa uncompressed.gfa --coverage simplified.wip.cov \
        --min-alt-ovl 10000 --max-unique-cov 35. -l 20000 -d 10000 --max-cov-ratio 1.5 --max-shortening 50 &> simple_unique_br_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Removing two-sided nongenomic links"
    $algo_root/nongenomic_link_removal simplified.wip.gfa uncompressed.gfa --coverage simplified.wip.cov \
        --unique-len 100000 --max-unique-cov 40. --reliable-cov 12. --reliable-len 20000 --both-sides &> two_sided_nongenomic_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Assigning coverage"
    $scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Removing one-sided nongenomic links"
    $algo_root/nongenomic_link_removal simplified.wip.gfa uncompressed.gfa --coverage simplified.wip.cov \
        --unique-len 200000 --max-unique-cov 40. --reliable-cov 12. --reliable-len 50000 &> one_sided_nongenomic_${final_it}.log
        #--unique-len 200000 --max-unique-cov 35. --reliable-cov 12. --reliable-len 50000 &> one_sided_nongenomic_${final_it}.log
    echo "Compressing round $cnt"
    $scripts_root/compact_gfa.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    cp simplified.wip.gfa simplified.final_${final_it}.gfa

done

cp simplified.wip.gfa simplified.gfa

$scripts_root/resolve_layouts.py simplified.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.txt
$scripts_root/assign_coverage.py resolved_mapping.txt $read_cov > simplified.cov
mv simplified.gfa no_cov.gfa
$scripts_root/inject_coverage.py no_cov.gfa simplified.cov > simplified.gfa
rm no_cov.gfa

rm -f *.wip.*

echo "All done"

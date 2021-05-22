#!/bin/bash
set -eou

if [ "$#" -lt 6 ]; then
    echo "script.sh <gfa> <utg mapping> <utg reads (.gfa with a-lines)> <read cov> <bubble_diff> <out_folder>"
    exit 239
fi

gfa=$(readlink -e $1)
utg_mapping=$(readlink -e $2)
utg_reads=$(readlink -e $3)
read_cov=$(readlink -e $4)
#bubble_diff=2000
bubble_diff=$5
out_folder=$6

scripts_root=$(dirname $(readlink -e $0))
#algo_root=$scripts_root/../gfacpp/build
algo_root=~/git/gfacpp/build/

mkdir -p $out_folder
cd $out_folder

echo "Processing graph $gfa"
echo "Bubble diff length set to $bubble_diff"

grep -v "^a" $gfa | grep -v "^x" > simplified.wip.gfa

cp $utg_mapping mapping.txt

#Compaction rounds counter
cnt=1

#Added to fix invalid overlaps
$algo_root/test simplified.wip.gfa simplified.wip.gfa --prefix f$((cnt++))_ --id-mapping mapping.txt &> normalize.log

#max_tip_len=25000
max_tip_len=20000
#max_lc_len=35000
max_lc_len=30000
unique_cov_thr=35

final_it_cnt=2

for final_it in $(seq 1 $final_it_cnt) ; do

    #weak_thr=$((4000 + final_it * 2000))
    weak_thr=$((4000 + final_it * 1000))
    echo "Removing overlaps weaker $weak_thr with dead-end prevention"
    $algo_root/weak_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --min-overlap $weak_thr --prevent-deadends &> weak_removal${final_it}.log

    echo "Removing tips"
    $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --max-length $max_tip_len &> tip_clipper${final_it}.log

    cov_thr=$((final_it + 3))

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    #echo "Assigning coverage"
    #$scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov > simplified.wip.cov
    echo "Removing low coverage nodes iteration $final_it (removing nodes with coverage < $cov_thr)"
    #$algo_root/low_cov_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage simplified.wip.cov --cov-thr $cov_thr --max-length 35000 &> low_cov${final_it}.log
    $algo_root/low_cov_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --cov-thr $cov_thr --max-length $max_lc_len &> low_cov${final_it}.log

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Killing loops"
    $algo_root/loop_killer simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --max-base-cov $unique_cov_thr &> loop_killer_${final_it}.log

    echo "Removing simple bulges ${final_it}"
    #$algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --max-length 10000 --max-diff 3 --min-alt-ovl 10000 &> simple_br_${final_it}.log
    $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --max-length 0 --max-diff 3 --min-alt-ovl 10000 &> simple_br_${final_it}.log

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Removing shortcuts"
    $algo_root/shortcut_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --max-base-cov $unique_cov_thr --min-path-cov 8 &> shortcut_removal_${final_it}.log

    echo "Final bubble removal ${final_it}"
    $algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --max-length 0 --max-diff $bubble_diff &> final_bfinder_${final_it}.log
    #$algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --max-length 60000 --max-diff $bubble_diff &> final_bfinder_${final_it}.log

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Removing low frequency hets in unique areas"
    $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) \
        --min-alt-ovl 10000 --max-unique-cov $unique_cov_thr --max-cov-ratio 0.33 --max-length 10000 --max-diff 500 &> low_freq_br_${final_it}.log
        #--min-alt-ovl 10000 --max-unique-cov 30. --max-cov-ratio 0.33 --max-length 30000 --max-diff 5000 &> low_freq_br_${final_it}.log

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Removing other variants in unique areas"
    $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) \
        --min-alt-ovl 10000 --max-unique-cov 28. --max-length 10000 --max-diff 500 --max-cov-ratio 1.5 --max-shortening 50 &> simple_unique_br_${final_it}.log
    #$algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) \
    #    --min-alt-ovl 10000 --max-unique-cov 35. --max-length 20000 --max-diff 10000 --max-cov-ratio 1.5 --max-shortening 50 &> simple_unique_br_${final_it}.log

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Removing two-sided nongenomic links"
    $algo_root/nongenomic_link_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) \
        --unique-len 150000 --max-unique-cov $unique_cov_thr --reliable-cov 15. --reliable-len 20000 --both-sides &> two_sided_nongenomic_${final_it}.log
    #$algo_root/nongenomic_link_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) \
    #    --unique-len 100000 --max-unique-cov 40. --reliable-cov 12. --reliable-len 20000 --both-sides &> two_sided_nongenomic_${final_it}.log

    #echo "Resolving layout"
    #$scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    #echo "Removing one-sided nongenomic links"
    #$algo_root/nongenomic_link_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix f$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) \
    #    --unique-len 200000 --max-unique-cov 40. --reliable-cov 12. --reliable-len 50000 &> one_sided_nongenomic_${final_it}.log
    #    #--unique-len 200000 --max-unique-cov 35. --reliable-cov 12. --reliable-len 50000 &> one_sided_nongenomic_${final_it}.log

    cp simplified.wip.gfa simplified.final_${final_it}.gfa

done

cp simplified.wip.gfa simplified.gfa
rm -f *.wip.*

$scripts_root/resolve_layouts.py simplified.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.txt
$scripts_root/assign_coverage.py resolved_mapping.txt $read_cov > simplified.cov
mv simplified.gfa no_cov.gfa
$scripts_root/inject_coverage.py no_cov.gfa simplified.cov > simplified.gfa
rm no_cov.gfa

echo -e "H\tVN:Z:1.0" > simplified.noseq.gfa
$scripts_root/../gfacpp/gfatools/gfatools view -S simplified.gfa >> simplified.noseq.gfa

$scripts_root/resolve_layouts.py simplified.noseq.gfa mapping.txt --partial-resolve > resolved_utg.txt

awk '/^S/{print ">"$2"\n"$3}' simplified.gfa | fold > simplified.nodes.fasta

echo "All done"
cd -

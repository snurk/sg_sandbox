#!/bin/bash
set -eou

if [ "$#" -lt 3 ]; then
    echo "script.sh <gfa> <utg reads (.gfa with a-lines)> <read cov>"
    exit 239
fi

gfa=$1
utg_reads=$2
read_cov=$3

scripts_root=$(dirname $(readlink -e $0))
algo_root=$scripts_root/../gfacpp/build

echo "Pruning graph $gfa"

if [ -n "$NO_DEADEND_OVERLAP_THRESHOLDS" ]; then
    echo "Final simplification rounds $NO_DEADEND_OVERLAP_THRESHOLDS"
else
    echo "Final simplification rounds have not been configured!!!"
    exit 239
fi

if [ -n "$WEAK_OVERLAP_THRESHOLDS" ]; then
    echo "Basic simplification rounds $WEAK_OVERLAP_THRESHOLDS"
else
    echo "Basic simplification rounds haven't been configured!!!"
    exit 239
fi

grep -v "^a" $gfa | grep -v "^x" > simplified.wip.gfa
rm -f mapping.txt
touch mapping.txt

#Compaction rounds counter
cnt=1

#Overlap normalization step
$algo_root/test simplified.wip.gfa simplified.wip.gfa --prefix m$((cnt++))_ --id-mapping mapping.txt &> normalize.log

for read_cnt_bound in $(seq 1 $S_MAX_INIT_TIP_READ_CNT) ; do

    echo "Read count aware tip clipping (read count bound $read_cnt_bound)"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> cnt_aware_tclipper_$read_cnt_bound.log
    $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --read-cnt-file <(sed 's/ /,/g' resolved_mapping.wip.txt | awk -F"," '{print $1,NF-1}') --max-read-cnt $read_cnt_bound --max-length $S_MAX_TIP_LEN &>> cnt_aware_tclipper_$read_cnt_bound.log

done

echo "Initial bubble removal"
# remove bubbles
$algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $S_MAX_BUBBLE_LEN --max-diff $S_BUBBLE_DIFF &> bfinder_init.log

weak_it=1
cp simplified.wip.gfa iteration0.gfa

echo "Will use iterative weak_ovl thresholds $WEAK_OVERLAP_THRESHOLDS"
for it_setting in $(echo $WEAK_OVERLAP_THRESHOLDS | tr ' ' '\n') ; do

    weak_ovl=$(echo $it_setting | awk -F: '{print $1}')
    cov_thr=$(echo $it_setting | awk -F: '{print $2}')
    if [ -z $cov_thr ] ; then
        echo "Wrong format of iteration thresholds configuration"
        exit 239
    fi

    echo "Removing weak connections with threshold $weak_ovl"
    # remove weak links if alternatives present
    $algo_root/weak_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --min-overlap $weak_ovl &> weak_removal${weak_it}.log

    echo "Clipping tips iteration $weak_it"
    $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $S_MAX_TIP_LEN &> tclipper${weak_it}.log

    echo "Removing bubbles iteration $weak_it"
    # remove bubbles
    $algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $S_MAX_BUBBLE_LEN --max-diff $S_BUBBLE_DIFF &> bfinder${weak_it}.log

    echo "Minimal coverage threshold set at $cov_thr"

    echo "Removing low coverage nodes iteration $weak_it (removing nodes with coverage < $cov_thr)"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> low_cov${weak_it}.log
    $algo_root/low_cov_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --cov-thr $cov_thr --max-length $S_MAX_LC_LEN &>> low_cov${weak_it}.log

    #echo "Removing links between nodes of widely different coverage $weak_it"
    #$scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> unbalanced${weak_it}.log
    #$algo_root/unbalanced_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --cov-ratio 0.1 &>> unbalanced${weak_it}.log

    cp simplified.wip.gfa iteration${weak_it}.${weak_ovl}.gfa
    weak_it=$((weak_it+1))
done
echo "Basic simplification done"

echo "Will use iterative no deadend thresholds $NO_DEADEND_OVERLAP_THRESHOLDS"
for it_setting in $(echo $NO_DEADEND_OVERLAP_THRESHOLDS | tr ' ' '\n') ; do

    weak_ovl=$(echo $it_setting | awk -F: '{print $1}')
    cov_thr=$(echo $it_setting | awk -F: '{print $2}')
    if [ -z $cov_thr ] ; then
        echo "Wrong format of iteration thresholds configuration"
        exit 239
    fi

    echo "Removing overlaps weaker $weak_ovl with dead-end prevention"
    $algo_root/weak_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --min-overlap $weak_ovl --prevent-deadends &> weak_removal${weak_it}.log

    #TODO introduce separate final tip length threshold?
    echo "Removing tips"
    $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $S_MAX_TIP_LEN &> final_tclipper${weak_it}.log

    echo "Removing low coverage nodes iteration $weak_it (removing nodes with coverage < $cov_thr)"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> low_cov${weak_it}.log
    $algo_root/low_cov_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --cov-thr $cov_thr --max-length $S_MAX_LC_LEN &>> low_cov${weak_it}.log

    echo "Killing loops"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> loop_killer_${weak_it}.log
    $algo_root/loop_killer simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --max-base-cov $S_UNIQUE_COV_THR &>> loop_killer_${weak_it}.log

    echo "Removing simple bulges ${weak_it}"
    $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $S_MAX_BUBBLE_LEN --max-diff $S_BUBBLE_DIFF --min-alt-ovl $S_MIN_ALT_OVL &> simple_br_${weak_it}.log

    echo "Removing shortcuts"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> shortcut_removal_${weak_it}.log
    $algo_root/shortcut_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --max-base-cov $S_UNIQUE_COV_THR --min-path-cov $S_RELIABLE_COV &>> shortcut_removal_${weak_it}.log

    echo "Final bubble removal ${weak_it}"
    $algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $S_MAX_BUBBLE_LEN --max-diff $S_BUBBLE_DIFF &> final_bfinder_${weak_it}.log

    if [ "${S_ADVANCED_FINAL_SIMPLIF}" = 0 ] ; then
        echo "Not using advanced final simplification procedures"
    else
        echo "Invoking advanced final simplification procedures"
        #FIXME parameterize
        echo "Removing low frequency hets in unique areas"
        $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> low_freq_br_${weak_it}.log
        $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --min-alt-ovl $S_MIN_ALT_OVL --max-unique-cov $S_UNIQUE_COV_THR --max-cov-ratio 0.33 --max-length 10000 --max-diff 500 &>> low_freq_br_${weak_it}.log

        #FIXME parameterize
        #FIXME isn't it still quite risky?
        echo "Removing other variants in unique areas"
        $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> simple_unique_br_${weak_it}.log
        $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --min-alt-ovl $S_MIN_ALT_OVL --max-unique-cov $S_CONS_UNIQUE_COV_THR --max-length 10000 --max-diff 500 --max-cov-ratio 1.5 --max-shortening 50 &>> simple_unique_br_${weak_it}.log

        #FIXME parameterize
        echo "Removing two-sided nongenomic links"
        $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> two_sided_nongenomic_${weak_it}.log
        $algo_root/nongenomic_link_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --unique-len $S_UNIQUE_LEN --max-unique-cov $S_UNIQUE_COV_THR --reliable-cov $S_RELIABLE_COV --reliable-len $S_RELIABLE_LEN --both-sides &>> two_sided_nongenomic_${weak_it}.log

        #echo "Removing one-sided nongenomic links"
        #$scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt 2> one_sided_nongenomic_${weak_it}.log
        #$algo_root/nongenomic_link_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --unique-len $S_UNIQUE_LEN --max-unique-cov $S_CONS_UNIQUE_COV_THR --reliable-cov $S_RELIABLE_COV --reliable-len $S_CONS_RELIABLE_LEN &>> one_sided_nongenomic_${weak_it}.log
    fi

    cp simplified.wip.gfa iteration${weak_it}.${weak_ovl}.gfa
    weak_it=$((weak_it+1))
done
echo "Final simplification done"

echo "Removing isolated nodes"
# remove isolated nodes
$algo_root/isolated_remover simplified.wip.gfa simplified.gfa --max-length $S_MAX_ISOLATED_LEN &> iremover.log

#rm -f iteration*.gfa
rm -f *.wip.*

echo "Simplification done"

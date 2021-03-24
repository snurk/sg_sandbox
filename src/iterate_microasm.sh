#!/bin/bash
set -eou

if [ "$#" -lt 4 ]; then
    echo "script.sh <gfa> <utg reads (.gfa with a-lines)> <read cov> <bubble_diff> [weak_overlaps...(default: disabled)]"
    exit 239
fi

gfa=$1
utg_reads=$2
#bubble_diff=2000
read_cov=$3
bubble_diff=$4

scripts_root=$(dirname $(readlink -e $0))

#TODO update repo and use local path
algo_root=$scripts_root/../gfacpp/build
#algo_root=~/git/gfacpp/build/

echo "Processing graph $gfa"
echo "Bubble diff length set to $bubble_diff"

if [ "$#" -lt 5 ]; then
    echo "Weak overlap removal disabled"
else
    echo "Will use iterative weak_ovl thresholds" "${@:5}"
fi

grep -v "^a" $gfa | grep -v "^x" > simplified.wip.gfa
rm -f mapping.txt
touch mapping.txt

#Compaction rounds counter
cnt=1

for read_cnt_bound in 1 2 3 ; do

    #FIXME crazy inefficient
    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    #To get read counts
    sed 's/ /,/g' resolved_mapping.wip.txt | awk -F"," '{print $1,NF-1}' > read_cnt.wip.txt
    echo "Read count aware tip clipping (read count bound $read_cnt_bound)"
    # clip tips
    #$algo_root/cnt_aware_tip_clipper simplified.wip.gfa uncompressed.gfa read_cnt.wip.txt $read_cnt_bound 25000 &> cnt_aware_tclipper_$read_cnt_bound.log
    $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --read-cnt-file read_cnt.wip.txt --max-read-cnt $read_cnt_bound --max-length 25000 &> cnt_aware_tclipper_$read_cnt_bound.log

done

echo "Initial bubble removal"
# remove bubbles
$algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length 60000 --max-diff $bubble_diff &> bfinder_init.log

weak_it=1
cp simplified.wip.gfa iteration0.gfa

if [ "$#" -lt 5 ]; then
    echo "Will not remove weak connections"
else
    echo "Will use iterative weak_ovl thresholds" "${@:5}"
    for weak_ovl in "${@:5}" ; do

        echo "Removing weak connections with threshold $weak_ovl"
        # remove weak links if alternatives present
        $algo_root/weak_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --min-overlap $weak_ovl &> weak_removal${weak_it}.log
        #$algo_root/weak_removal simplified.wip.gfa uncompressed.gfa $weak_ovl &> weak_removal${weak_it}.log

        if [ $weak_it -eq 1 ] ; then
            for init_tip_bound in 5000 15000 25000 ; do

                echo "Initial tip clipping (length bound $init_tip_bound)"
                # clip tips
                $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length $init_tip_bound &> tclipper_init_$init_tip_bound.log
                #$algo_root/tip_clipper simplified.wip.gfa uncompressed.gfa $init_tip_bound &> tclipper_init_$init_tip_bound.log

            done
        fi

        echo "Clipping tips iteration $weak_it"
        $algo_root/tip_clipper simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length 25000 &> tclipper${weak_it}.log
        # clip tips
        #$algo_root/tip_clipper simplified.wip.gfa uncompressed.gfa 25000 &> tclipper${weak_it}.log

        echo "Removing bubbles iteration $weak_it"
        # remove bubbles
        $algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length 60000 --max-diff $bubble_diff &> bfinder${weak_it}.log
        #$algo_root/bubble_removal simplified.wip.gfa uncompressed.gfa 60000 $bubble_diff &> bfinder${weak_it}.log

        if [ $weak_ovl -gt 7000 ] ; then
            cov_thr=5
        elif [ $weak_ovl -gt 5000 ] ; then
            cov_thr=4
        elif [ $weak_ovl -gt 3000 ] ; then
            cov_thr=3
        else
            cov_thr=2
        fi

        echo "Minimal coverage threshold set at $cov_thr"

        echo "Resolving layout"
        $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
        echo "Removing low coverage nodes iteration $weak_it (removing nodes with coverage < $cov_thr)"
        #Remove everything of coverage < 4
        #Used to be 50000
        $algo_root/low_cov_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --cov-thr $cov_thr --max-length 30000 &> low_cov${weak_it}.log
        #$algo_root/low_cov_remover simplified.wip.gfa uncompressed.gfa simplified.wip.cov $cov_thr 30000 &> low_cov${weak_it}.log

        #echo "Resolving layout"
        #$scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
        #echo "Removing links between nodes of widely different coverage $weak_it"
        #$algo_root/unbalanced_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --cov-ratio 0.1 &> unbalanced${weak_it}.log

        cp simplified.wip.gfa iteration${weak_it}.gfa

        ln -sf iteration${weak_it}.gfa simplified.$weak_ovl.gfa

        weak_it=$((weak_it+1))
    done
    echo "Iterative weak removal done"
fi

cp mapping.txt pre_final.mapping.txt
cp simplified.wip.gfa pre_final.simplified.gfa

final_it_cnt=3

for final_it in $(seq 1 $final_it_cnt) ; do

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Killing loops"
    $algo_root/loop_killer simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --max-base-cov 30 &> loop_killer_${final_it}.log
    #$algo_root/loop_killer simplified.wip.gfa uncompressed.gfa simplified.wip.cov 30 &> loop_killer_${final_it}.log

    echo "Removing simple bulges ${final_it}"
    $algo_root/simple_bulge_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length 10000 --max-diff 3 --min-alt-ovl 10000 &> simple_br_${final_it}.log
    #$algo_root/simple_bulge_removal simplified.wip.gfa uncompressed.gfa -l 10000 -d 3 --min-alt-ovl 10000 &> simple_br_${final_it}.log

    echo "Resolving layout"
    $scripts_root/resolve_layouts.py simplified.wip.gfa mapping.txt --miniasm $utg_reads > resolved_mapping.wip.txt
    echo "Removing shortcuts"
    $algo_root/shortcut_remover simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --coverage <($scripts_root/assign_coverage.py resolved_mapping.wip.txt $read_cov) --max-base-cov 30 --min-path-cov 8 &> shortcut_removal_${final_it}.log
    #$algo_root/shortcut_remover simplified.wip.gfa uncompressed.gfa simplified.wip.cov 30 8 &> shortcut_removal_${final_it}.log

    echo "Final bubble removal ${final_it}"
    $algo_root/bubble_removal simplified.wip.gfa simplified.wip.gfa --compact --prefix m$((cnt++))_ --id-mapping mapping.txt --max-length 60000 --max-diff $bubble_diff &> final_bfinder_${final_it}.log
    #$algo_root/bubble_removal simplified.wip.gfa uncompressed.gfa 60000 $bubble_diff &> final_bfinder_${final_it}.log

    cp simplified.wip.gfa simplified.final_${final_it}.gfa

done

echo "Removing isolated nodes"
# remove isolated nodes
$algo_root/isolated_remover simplified.wip.gfa simplified.gfa --max-length 20000 &> iremover.log

#rm -f iteration*.gfa
rm -f *.wip.*

echo "All done"

#!/bin/bash
set -e

if [ "$#" -lt 2 ]; then
    echo "script.sh <gfa> <bubble_diff> [weak_overlaps...(default: disabled)]"
    exit 239
fi

gfa=$1
#bubble_diff=2000
bubble_diff=$2

echo "Processing graph $gfa"
echo "Bubble diff length set to $bubble_diff"

if [ "$#" -lt 3 ]; then
    echo "Weak overlap removal disabled"
else
    echo "Will use iterative weak_ovl thresholds" "${@:3}"
fi

grep -v "^a" $gfa | grep -v "^x" > simplified.wip.gfa
rm -f mapping.txt
touch mapping.txt

#Compaction rounds counter
cnt=1

for read_cnt_bound in 1 2 3 ; do

    #FIXME crazy inefficient
    echo "Resolving layout"
    ~/git/ngs_scripts/bogart_gfa/resolve_mapping.py simplified.wip.gfa microasm.gfa mapping.txt > resolved_mapping.wip.txt
    #To get read counts
    sed 's/ /,/g' resolved_mapping.wip.txt | awk -F"," '{print $1,NF-1}' > read_cnt.wip.txt
    echo "Read count aware tip clipping (read count bound $read_cnt_bound)"
    # clip tips
    ~/git/ngs_scripts/gfakluge/cnt_aware_tip_clipper simplified.wip.gfa uncompressed.gfa read_cnt.wip.txt $read_cnt_bound 25000 &> cnt_aware_tclipper_$read_cnt_bound.log
    echo "Compressing round $cnt"
    # compress
    #Doesn't work for isolated segments
    #~/git/gfatools/gfatools asm -u uncompressed.gfa > simplified.gfa
    ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

done

echo "Initial bubble removal"
# remove bubbles
~/git/gfacpp/bubble_removal simplified.wip.gfa uncompressed.gfa 60000 $bubble_diff &> bfinder_init.log
echo "Compressing round $cnt"
~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
cnt=$((cnt+1))
rm uncompressed.gfa

weak_it=1
cp simplified.wip.gfa iteration0.gfa

if [ "$#" -lt 3 ]; then
    echo "Will not remove weak connections"
else
    echo "Will use iterative weak_ovl thresholds" "${@:3}"
    for weak_ovl in "${@:3}" ; do

        echo "Removing weak connections with threshold $weak_ovl"
        # remove weak links if alternatives present
        ~/git/gfacpp/weak_removal simplified.wip.gfa uncompressed.gfa $weak_ovl &> weak_removal${weak_it}.log
        # FIXME get rid of L1/L2 output
        #sed -i 's/L1/l1/g' uncompressed.gfa
        #sed -i 's/L2/l2/g' uncompressed.gfa
        #cp uncompressed.gfa no_weak${weak_it}.nocompress.gfa
        echo "Compressing round $cnt"
        ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
        cnt=$((cnt+1))
        rm uncompressed.gfa

        if [ $weak_it -eq 1 ] ; then
            for init_tip_bound in 5000 15000 25000 ; do

                echo "Initial tip clipping (length bound $init_tip_bound)"
                # clip tips
                ~/git/ngs_scripts/gfakluge/tip_clipper simplified.wip.gfa uncompressed.gfa $init_tip_bound &> tclipper_init_$init_tip_bound.log
                echo "Compressing round $cnt"
                # compress
                #Doesn't work for isolated segments
                #~/git/gfatools/gfatools asm -u uncompressed.gfa > simplified.gfa
                ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
                cnt=$((cnt+1))
                rm uncompressed.gfa

            done
        fi

        echo "Clipping tips iteration $weak_it"
        # clip tips
        ~/git/ngs_scripts/gfakluge/tip_clipper simplified.wip.gfa uncompressed.gfa 25000 &> tclipper${weak_it}.log
        echo "Compressing round $cnt"
        ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
        cnt=$((cnt+1))
        rm uncompressed.gfa

        echo "Removing bubbles iteration $weak_it"
        # remove bubbles
        ~/git/gfacpp/bubble_removal simplified.wip.gfa uncompressed.gfa 60000 $bubble_diff &> bfinder${weak_it}.log
        echo "Compressing round $cnt"
        ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
        cnt=$((cnt+1))
        rm uncompressed.gfa

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
        ~/git/ngs_scripts/bogart_gfa/resolve_mapping.py simplified.wip.gfa microasm.gfa mapping.txt > resolved_mapping.wip.txt
        echo "Assigning coverage"
        ~/git/ngs_scripts/gfakluge/assign_coverage.py resolved_mapping.wip.txt min_read.cov > simplified.wip.cov
        echo "Removing low coverage nodes iteration $weak_it (removing nodes with coverage < $cov_thr)"
        #Remove everything of coverage < 4
        #Used to be 50000
        ~/git/ngs_scripts/gfakluge/low_cov_remover simplified.wip.gfa uncompressed.gfa simplified.wip.cov $cov_thr 30000 &> low_cov${weak_it}.log
        echo "Compressing round $cnt"
        ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
        cnt=$((cnt+1))
        rm uncompressed.gfa

        #echo "Resolving layout"
        #~/git/ngs_scripts/bogart_gfa/resolve_mapping.py simplified.wip.gfa microasm.gfa mapping.txt > resolved_mapping.wip.txt
        #echo "Assigning coverage"
        #~/git/ngs_scripts/gfakluge/assign_coverage.py resolved_mapping.wip.txt min_read.cov > simplified.wip.cov
        #echo "Removing links between nodes of widely different coverage $weak_it"
        #~/git/gfacpp/unbalanced_removal simplified.wip.gfa uncompressed.gfa simplified.wip.cov 0.1 &> unbalanced${weak_it}.log
        #echo "Compressing round $cnt"
        #~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
        #cnt=$((cnt+1))
        #rm uncompressed.gfa

        cp simplified.wip.gfa iteration${weak_it}.gfa

        ln -s iteration${weak_it}.gfa simplified.$weak_ovl.gfa

        weak_it=$((weak_it+1))
    done
    echo "Iterative weak removal done"
fi

final_it_cnt=3

for final_it in $(seq 1 $final_it_cnt) ; do

    echo "Resolving layout"
    ~/git/ngs_scripts/bogart_gfa/resolve_mapping.py simplified.wip.gfa microasm.gfa mapping.txt > resolved_mapping.wip.txt
    echo "Assigning coverage"
    ~/git/ngs_scripts/gfakluge/assign_coverage.py resolved_mapping.wip.txt min_read.cov > simplified.wip.cov
    echo "Killing loops"
    ~/git/gfacpp/loop_killer simplified.wip.gfa uncompressed.gfa simplified.wip.cov 30 > loop_killer_${final_it}.log
    echo "Compressing round $cnt"
    ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Removing simple bulges ${final_it}"
    ~/git/gfacpp/simple_bulge_removal simplified.wip.gfa uncompressed.gfa 10000 3 10000 &> simple_br_${final_it}.log
    echo "Compressing round $cnt"
    ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Resolving layout"
    ~/git/ngs_scripts/bogart_gfa/resolve_mapping.py simplified.wip.gfa microasm.gfa mapping.txt > resolved_mapping.wip.txt
    echo "Assigning coverage"
    ~/git/ngs_scripts/gfakluge/assign_coverage.py resolved_mapping.wip.txt min_read.cov > simplified.wip.cov
    echo "Removing shortcuts"
    ~/git/gfacpp/shortcut_remover simplified.wip.gfa uncompressed.gfa simplified.wip.cov 30 8 &> shortcut_removal_${final_it}.log
    ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    echo "Final bubble removal ${final_it}"
    ~/git/gfacpp/bubble_removal simplified.wip.gfa uncompressed.gfa 60000 $bubble_diff &> final_bfinder_${final_it}.log
    echo "Compressing round $cnt"
    ~/git/ngs_scripts/bogart_gfa/compress.py uncompressed.gfa simplified.wip.gfa m${cnt}_ 2>> mapping.txt
    cnt=$((cnt+1))
    rm uncompressed.gfa

    cp simplified.wip.gfa simplified.final_${final_it}.gfa

done

echo "Removing isolated nodes"
# remove isolated nodes
~/git/ngs_scripts/gfakluge/isolated_remover simplified.wip.gfa simplified.gfa &> iremover.log

#rm -f iteration*.gfa
rm -f *.wip.*

echo "All done"
#TODO use gimbricate before microasm!?
#~/git/gimbricate/bin/gimbricate -g simplified.gfa -n -t 8 > refined_ovl.gfa

#For bluntification?
#seqwish -s h.fasta -p h.paf -g h.seqwish.gfa

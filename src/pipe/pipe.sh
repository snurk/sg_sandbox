#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "script.sh <canu bin> <canu assembly folder> [min_overlap = 2000] [max_erate = 0.001]"
    echo "NB: for HiFi analysis reads should be in homopolymer compressed space to match the overlaps"
    echo "NB2!!!: never use max_erate higher than what was used by bogart"
    exit 239
fi

scripts_root=$(dirname $(readlink -e $0))
bin=$1
assembly=$2

if [ "$#" -lt 3 ]; then
    min_ovl=2000
else
    min_ovl=$3
    echo "Will use min_ovl threshold of $min_ovl"
fi

if [ "$#" -lt 4 ]; then
    max_erate=0.001
else
    max_erate=$4
    echo "Will use overlap error threshold of $max_erate"
fi

if [ -f reads.compressed.fasta.gz ]; then
    echo "File reads.compressed.fasta.gz was found in the folder and will be reused"
elif [ -f ../reads.compressed.fasta.gz ]; then
    echo "File reads.compressed.fasta.gz was found in the folder and will be linked"
    ln -sf ../reads.compressed.fasta.gz
else
    #-compressed -trimmed
    echo "Getting reads from assembly folder"
    $bin/sqStoreDumpFASTQ -S $assembly/asm.seqStore -nolibname -noreadname -fasta -o reads.compressed.wip.gz
    mv reads.compressed.wip.fasta.gz reads.compressed.fasta.gz
fi

ovlstore=$assembly/unitigging/4-unitigger/asm.0.all.ovlStore
#ovlstore=$assembly/unitigging/asm.ovlStore

if [ -f min_read.cov ]; then
    echo "File min_read.cov was found in the folder and will be reused"
elif [ -f ../min_read.cov ]; then
    echo "File min_read.cov was found in the root folder and will be linked"
    ln -sf ../min_read.cov
else
    echo "Getting read coverage estimates"
    $bin/ovStoreDump -S $assembly/asm.seqStore -O $ovlstore -coverage -erate 0-$max_erate | tail -n+4 | awk '{print $1,$2}' > min_read.wip.cov
    sed '/^\s*$/d' min_read.wip.cov > min_read.cov
    rm -f min_read.wip.cov
fi

if [ -f simplified.gfa ]; then
    echo "File simplified.gfa was found in the folder and will be reused"
else

    if [ -f ovl.paf ]; then
        echo "File ovl.paf was found in the folder and will be reused"
    else
        echo "Getting overlaps"
        #-nobogartspur
        $bin/ovStoreDump -bogart $assembly/unitigging/4-unitigger/asm.best.edges -nobogartcontained -nobogartcoveragegap -noredundant -nocontained -nocontainer -erate 0-$max_erate -paf -S $assembly/asm.seqStore/ -O $ovlstore | awk 'BEGIN { OFS = "\t"} {$1="read"$1; $6="read"$6; print $0}' > ovl.tmp.paf && mv ovl.tmp.paf ovl.paf
    fi

    if [ -f processed.gfa ]; then
        echo "File processed.gfa was found in the folder and will be reused"

        if [ ! -f microasm.gfa ]; then
            echo "File microasm.gfa is missing"
            exit 239
        fi
    else
        echo "Building initial graph"
        $scripts_root/build_graph.sh reads.compressed.fasta.gz ovl.paf $min_ovl
    fi

    grep "^a" microasm.gfa > utg_reads.gfa
    $scripts_root/simplif.sh processed.gfa utg_reads.gfa min_read.cov
fi

if [ ! -f microasm.gfa ]; then
    echo "File microasm.gfa is missing"
    exit 239
fi

if [ ! -f mapping.txt ]; then
    echo "File mapping.txt is missing"
    exit 239
fi

awk '/^S/{print ">"$2"\n"$3}' simplified.gfa | fold > simplified.nodes.fasta

$scripts_root/../resolve_layouts.py simplified.gfa mapping.txt --miniasm utg_reads.gfa > resolved_mapping.txt

$scripts_root/../assign_coverage.py resolved_mapping.txt min_read.cov > simplified.cov

mv simplified.gfa no_cov.gfa
$scripts_root/../inject_coverage.py simplified.cov < no_cov.gfa > simplified.gfa
rm no_cov.gfa

echo -e "H\tVN:Z:1.0" > simplified.noseq.gfa
$scripts_root/../../gfacpp/gfatools/gfatools view -S simplified.gfa >> simplified.noseq.gfa

$scripts_root/../resolve_layouts.py simplified.noseq.gfa mapping.txt --partial-resolve > resolved_utg.txt

rm -f utg_reads.gfa

echo "Pipeline done"

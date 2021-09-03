#!/bin/bash

if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <gap name>"
    exit 1
fi


gap=$1

echo "***** Processing $1 ******"
for i in `ls  ${gap}v*fasta`; do
   seqtk comp $i |awk '{print $1"\t"$2}'
done

max=0
best=()
for i in `ls $1v*.paf`; do
   good=`cat $i |awk '{if ($4 > 90 && $6 < 500 && $7 +500 > $8) print $0}'|wc -l`
   if [ $good -gt $max ]; then
      max=$good
      best=("$i")
   elif [ $good -eq $max ]; then
      best+=("$i")
   fi

   echo "**** Checking version $i and it has $good vs $max best"
   cat $i |grep "#"
   cat $i |awk '{if ($4 > 90 && $3 < -5000 && ($6 > 500 || $7 +500 < $8) && $10 > 500 && $11+500  < $12) print $0}' 
done

if [ ${#best[@]} -eq 1 ]; then
   echo "Best score for gap $gap was $max, winner is ${best[0]}"
else
   echo "Best score for gap $gap was $max, ${#best[@]} winners, comparing identities"
   e=$((${#best[@]}-1))

   for i in `seq 0 $e`; do
      s=$(($i+1))
      for j in `seq $s $e`; do
         echo "**** Comparing ${best[$i]} vs ${best[$j]}"
         cat ${best[$i]} |awk '{if ($4 > 90 && $6 < 500 && $7 +500 > $8) print $1"\t"$4}'|sort -k1 > tmp
         cat ${best[$j]} |awk '{if ($4 > 90 && $6 < 500 && $7 +500 > $8) print $1"\t"$4}'|sort -k1 > tmp2
         join tmp tmp2 |awk -v F="${best[$i]}" -v S="${best[$j]}" '{if ($2 > $3) print F;  else if ($3 > $2) print S; else print "TIE"; }' |sort |uniq -c
       done
   done
fi

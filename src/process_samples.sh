#!/bin/bash

# CALL: bash process_samples.sh /media/ultron/2tb_disk2/0_startallover/followup_meQTLs/cord_blood/test2/

DIR=$1
cd $DIR
mkdir PROCESSED
FILES=($( ls *.txt -1v))
l=${#FILES[@]}

for ((i=0;i<=$l-1;i+=1))
do
 # Combine Cs in Fwd/Rv to make up a CpG site
 index="$(cut -d'.' -f1 <<< ${FILES[$i]})" 
 echo ${FILES[$i]}
 echo $index
 #paste -d " "  - - < ${FILES[$i]} | awk '($4 > 0) || ($5 > 0) || ($11 > 0) || ($12 > 0)' > ./PROCESSED/${index}_processed.txt
 awk -F'\t' '$3 == "+"' ${FILES[$i]} | awk -F'\t' '$3=$2+1' | sed 's/ /\t/g' > tmp_p.bed
 awk -F'\t' '$3 == "-"' ${FILES[$i]} | awk -F'\t' '$3=$2, $2=($2-1)'| sed 's/ /\t/g' > tmp_m.bed
 bedtools intersect -a tmp_p.bed -b tmp_m.bed -wa -wb | awk '($4 > 0) || ($5 > 0) || ($11 > 0) || ($12 > 0)' > ./PROCESSED/${index}_processed.txt

 # Split file per CHR
 cd ./PROCESSED
 mkdir CHR_${index}
 mv ${index}_processed.txt ./CHR_${index}
 cd CHR_${index}
 awk -F '\t' '{print>$1}' ${index}_processed.txt
 rm ${index}_processed.txt
 cd ../../
done




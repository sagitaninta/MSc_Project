#!/bin/bash

## Change folder to run script in:
CHR=chr${1:48:-8} # Extract chromosome name

mkdir $CHR
mv $1 $CHR
cd $CHR

## Filter out cds and intersect with ancestral alleles:

gunzip -c $1 | awk '$3~"CDS"' | grep -v ncRNA_gene | grep -v pseudo | \
sort -k1,1 -k4,4n | bedtools merge -i - | bedops --chop - | sed 's/^/chr/' | \
bedtools getfasta -fi ../../../Pipeline_files/07_07_tests/AndeanFox_canFam3.1.fa -bed - -bedOut | \
grep [ATGC] > cds.bed

## Get sift scores for each alternative allele:
for i in A C G T;
do \
# Convert cds positions to vep format file:
sed 's/chr//' cds.bed | awk '{print $1 "\t" $2+1 "\t" $3 "\t" $4"/"}' | sed "s/$/$i/"  > ${i}_cds.bed

# Run vep:
time vep -i ${i}_cds.bed --offline --cache --dir ../../../miniconda3/bin/modules/Bio/EnsEMBL/VEP \
--species "canis_lupus_familiaris" --force_overwrite --sift b --tab -o ${i}_vep.txt \
--fields "Location,Allele,Consequence,SIFT"

# Extract high confidence sift scores from vep output file:
grep -v "#" ${i}_vep.txt | grep -v "-" | sed 's/:/\t/' | sed 's/(/\t/' | sed 's/)//' | \
grep -v "low_confidence" | grep -v "tolerated" | \
awk '{print $1 "\t" $2-1 "\t" $2 "\t" $3  "\t" $6}' > ${i}_sift_high.bed

# Extract lowest sift score (most deleterious) for each position:
awk '{print $1 "\t" $2 "\t" $3}' ${i}_sift_high.bed | uniq | \
bedtools map -a - -b ${i}_sift_high.bed -o absmin > ${i}_sift_min.bed;
done

## Combine sift scores into one file:

# Create bed file with all sift positions:
bedops -u *min.bed | awk '{print $1 "\t" $2 "\t" $3}' | uniq > all_sift_positions.bed

# Overlap sift scores for each allele with all positions (replace coloumns with no score with 1
# - will be 0 in calculation:

sed 's/^chr//' cds.bed | bedtools intersect -a all_sift_positions.bed -b - -sorted -wb | cut -f 1-3,7 | \
bedtools intersect -a - -b A_sift_min.bed -loj -sorted | cut -f 1-4,8 | \
bedtools intersect -a - -b C_sift_min.bed -loj -sorted | cut -f 1-5,9 | \
bedtools intersect -a - -b G_sift_min.bed -loj -sorted | cut -f 1-6,10 | \
bedtools intersect -a - -b T_sift_min.bed -loj -sorted | cut -f 1-7,11 | \
sed 's/\t\./\t1/g' > ${CHR}_sift_scores.bed


# Remove tmp files:
#rm cds.bed
#rm A_*
#rm C_*
#rm G_*
#rm T_*


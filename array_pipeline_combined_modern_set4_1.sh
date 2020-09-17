#!/bin/bash

#$ -cwd           #Run in current working directory  
#$ -V             #Make verbose
#$ -j y           #Join stdout and stderr
#$ -pe smp 1      #Request 1 core
#$ -l h_rt=24:0:0 #24 hour runtime
#$ -l h_vmem=3G   #3G per core (3GB RAM total)
#$ -t 1-10        #Number of bams/jobs to run
#$ -tc 3          # Task concurrency

module load angsd
module load python
module load bedtools
module load samtools

cd /data/scratch/bt19616/analysis_run_one/set4_modern/

BAM=$(sed -n "${SGE_TASK_ID}p" bam_list_set4_1.txt)

OUT=${BAM%.realigned_reduced.bam} # Extract sample name
SITES=/data/scratch/bt19616/analysis_run_one/ref_files/dog_angsd.file  # Specify co-ordinates for sites for all scores

## Specify score files:
AUTOCONS=/data/scratch/bt19616/analysis_run_one/ref_files/autosomes_anc_phylop_phastcons.bed
XCONS=/data/scratch/bt19616/analysis_run_one/ref_files/chrX_anc_phylop_phastcons.bed
AUTOSIFT=/data/scratch/bt19616/analysis_run_one/ref_files/autosomes_sift_scores.bed 
XSIFT=/data/scratch/bt19616/analysis_run_one/ref_files/chrX_sift_scores.bed

## Get genotype likelihoods:
samtools index $BAM

time angsd -i $BAM -sites $SITES -out $OUT -minQ 20 -minMapQ 20 -remove_bads 1 -trim 5 -GL 2 -doGlf 4

## Get genotype probabilities:

gunzip ${OUT}.glf.gz

source ~/python_env/bin/activate

time python ~/Masters_project/pipeline_scripts/all_genotype_likelihoods_v3.py ${OUT}.glf ${OUT}.gpf

rm ${OUT}.glf

## Get mutational load scores:
# Intersect probabilities with scores, get scores and remove bed file after use:
# Autosome cons scores:
bedtools intersect -a ${OUT}.gpf -b $AUTOCONS -sorted -wb | awk '{print $14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'  > ${OUT}_autocons.bed
time python ~/Masters_project/pipeline_scripts/mut_load_calculator_v3.py ${OUT}_autocons.bed ${OUT}_autocons
rm ${OUT}_autocons.bed

# Chr X cons scores:
bedtools intersect -a ${OUT}.gpf -b $XCONS  -sorted -wb | awk '{print $14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'  > ${OUT}_xcons.bed
time python ~/Masters_project/pipeline_scripts/mut_load_calculator_v3.py ${OUT}_xcons.bed ${OUT}_xcons
rm ${OUT}_xcons.bed

# Autosome sift scores:
bedtools intersect -a ${OUT}.gpf -b $AUTOSIFT  -sorted -wb | awk '{print $14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'  > ${OUT}_autosift.bed
time python ~/Masters_project/vep/sift_calculator_v2.py ${OUT}_autosift.bed ${OUT}_autosift_scores.txt
rm ${OUT}_autosift.bed

# Chr X sift scores:
bedtools intersect -a ${OUT}.gpf -b $XSIFT  -sorted -wb | awk '{print $14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}'  > ${OUT}_xsift.bed
time python ~/Masters_project/vep/sift_calculator_v2.py ${OUT}_xsift.bed ${OUT}_xsift_scores.txt
rm ${OUT}_xsift.bed

# Remove gpf file:
rm ${OUT}.gpf

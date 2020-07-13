import time
import sys
import csv
import numpy as np
import regex as re

#### Code ####

## Specify input file and extract sample name:
bed_file = sys.argv[1]
sample_name = re.split("\.", bed_file) # Split on dots in name
sample_name = sample_name[0] # Take name before first dot to get sample name

## Create dictionary of positions of homozygous transversions:

ATV = (4, 9) # CC and TT
CTV = (0, 7) # AA and GG
GTV = (4, 9) # CC and TT
TTV = (0, 7) # AA and GG

transversions = {"A":ATV, "C":CTV, "G":GTV, "T":TTV}

## Create dictionary of position of each homozygous genotype:

genotypes = {"A": 0, "C": 4, "G": 7, "T":9}

## Set total load scores and probability of homozygous positions to 0:

# Homozygous positions:
total_hom_pr = 0.0 # All homozygous positions
phylop_hom_pr = 0.0 # Homozygous positions with phylop score
phastcons_hom_pr = 0.0 # Homozygous positions with phastcons score 

# Load scores:
total_phylop = 0.0
total_phastcons = 0.0

## Calculate load for each line and pr of having homozygous derived/homozygous transversion:

with open(bed_file, "r") as file:
	bed = csv.reader(file, delimiter = "\t")
	for line in bed:
		if line[3] in ["A", "C", "G", "T"]:
		# Split line into scores, ancestral allele and posterior probabilities:
			post_prs = [float(i) for i in line[6:16]] # Split out posterior probabilities
			phylop = float(line[4]) # Extract phylop score
			phastcons = float(line[5]) # Extract phastcons score
			anc = line[3] # Extract ancestral genotype

		# Calculate load for homozygous transversions:
			pr_transv = 0.0
			for genotype in transversions[anc]:
				pr_transv += post_prs[genotype] # Add probability of homozygous transversions
			phylop_load = phylop * pr_transv
			phastcons_load = phastcons * pr_transv

		# Calculate probability  of position being ancestral homozygous or tranversion homozygous:
			hom_anc = genotypes[anc] # Co-ordinate of ancestral genotype in post_pr
			pr_anc =  post_prs[hom_anc] # Probability of ancestral genotype
			pr_hom = pr_transv + pr_anc # Probability transv or anc homozygous

		# Add load scores and probability of being homozygous to total:
			total_phylop += phylop_load
			total_phastcons += phastcons_load
			if phylop != 0:
				phylop_hom_pr += pr_hom  # Add pr if has a phylop score for position
			if phastcons != 0:
				phastcons_hom_pr += pr_hom # Add pr if has a phastcons score for position
			total_hom_pr += pr_hom # Add pr hom to hom total

## Write to file:
out_file = sample_name + "_load_scores.txt"

with open(out_file, "w") as file:
	(file.write(sample_name + "\t" + str(total_phylop) + "\t" + str(total_phastcons) + "\t" + 
	str(total_hom_pr) + "\t" + str(phylop_hom_pr) + "\t" + str(phastcons_hom_pr) + "\t" + 
	str(total_phylop/phylop_hom_pr) + "\t" + str(total_phastcons/phastcons_hom_pr) + "\n"))

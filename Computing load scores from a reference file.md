# Calculating load from pre-computed conservation scores

## Computing load scores from a reference file and posterior probability file

 Now we have posterior probabilties for our sample and a reference file containing conservation scores, we can compute
 a load score for our sample.  
 
### Step 1: Intersect posterior probability file with reference file  
 
Firstly we need to intersect the .gpf file and reference file so we have the posterior probabilities, phyloP scores, and 
phastcons scores in the same bed file:
```linux
bedtools intersect -a $SCORES -b ${OUT}.gpf -sorted -wb | cut -f 1-6,10-19 > ${OUT}.bed
```
This creates a bed file looking like this:
(add in image)

### Step 2: Calculate load scores using python

We can now use a python script that takes the bed file created in the last step and computes a load score across the
genome from this. The script essentially works by multiplying the probability of having a homozygous tranversion by
the conservation score at each position. The script also calculates the probability of observing a reference homozygous 
genotype for that position. These scores are then added to running totals. A break down of the script and the command
line argument used to run it.

1. Import modules:
```python
import time
import sys
import csv
import numpy as np
import regex as re
```

2. Specify input and outfiles from command line and extract sample name:
```python
bed_file = sys.argv[1]
out_file = sys.argv[2]
sample_name = re.split("\.", bed_file) # Split on dots in name
sample_name = sample_name[0] # Take name before first dot to get sample name
```

3. Set up dictionaries to extract genotypes:
a.) Dictionary of homozygous transversions
```python
ATV = (4, 9) # CC and TT
CTV = (0, 7) # AA and GG
GTV = (4, 9) # CC and TT
TTV = (0, 7) # AA and GG

transversions = {"A":ATV, "C":CTV, "G":GTV, "T":TTV}
```

b.) Dictionary of homozygous genotypes
genotypes = {"A": 0, "C": 4, "G": 7, "T":9}

4. Define total load scores and probability of observing homozygous reference/transversion and set to 0:
```python
total_hom_pr = 0.0 # All homozygous positions
phylop_hom_pr = 0.0 # Homozygous positions with phylop score
phastcons_hom_pr = 0.0 # Homozygous positions with phastcons score

total_phylop = 0.0
total_phastcons = 0.0
```

5. Calculate scores and probabilities for each line in the bed file:
```python
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
```
6. Write scores and number of homozygous positions to file:
```python
with open(out_file, "w") as file:
        (file.write(sample_name + "\t" + str(total_phylop) + "\t" + str(total_phastcons) + "\t" +
        str(total_hom_pr) + "\t" + str(phylop_hom_pr) + "\t" + str(phastcons_hom_pr) + "\t" +
        str(total_phylop/phylop_hom_pr) + "\t" + str(total_phastcons/phastcons_hom_pr) + "\n"))
```

7. Run the script from the command line:
```linux
time python 06_07_mut_load_calculator_v2.py ${OUT}_auto.bed ${OUT}_auto_load_scores.txt
```

### Step 3: Check the load score file:

The load score file generated above should look something like this:

If you have multiple samples you want to compare, these files can be combined using cat:



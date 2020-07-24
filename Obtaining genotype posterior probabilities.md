# Obtaining genotype posterior probabilities with ANGSD

Working with ancient DNA and low coverage samples can lead to uncertainty in genotype calls. One way to account for  
this uncertainty in the data is to maintain it in downstream analses and to weight any calculations by the genotype posterior  
probabilities, rather than doing a traditional genotype call. This can be done using ANGSD to obtain genotype likelihoods and  
then transforming these into posterior probabilities.

## Running ANGSD to get genotype likelihoods

### Step 1: Creating a sites file for ANGSD (if running on a subset of the genome)

ANGSD allows filtering by both chromsome (requires an indexed bam) or by particular sites. This requires an angsd file which  
can easily be created from a bed file using awk and then indexed with angsd http://www.popgen.dk/angsd/index.php/Sites 

1. Use awk to convert bed file with co-ordinates into an angsd site file:
```linux
awk '{print $1 "\t" $2+1 "\t" $3}' input.bed > angsd.file
```
2. Index the file:
```linux
angsd sites index angsd.file
```

### Step 2: Running ANGSD

ANGSD has various options for obtaining genotype likelihoods, all explained on the software's [wiki page](http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods). Below is 
an example command using a sites files and specifying a region (remember to index the bam first if using the -r option). This  
example is using a GATK model (-GL 2) and outputting all 10 genotype likelihoods in a tab separated format (-doGlf 4).

1. Specify files and run ANGSD with parameters of choice:
```linux
BAM=test_dog.bam
SITES=dog_angsd.file
REF=Canis_familiaris.CanFam3.1.dna.toplevel.fa 

angsd -i $BAM -sites $SITES -r chr10 -out $OUT -minQ 20 -minMapQ 20 -remove_bads 1 -ref $REF -GL 2 -doGlf 4
```

## Transforming genotype likelihoods into genotype probabilities

The step above will create a .glf file with a line for every position, with genotype likelihoods for each of the ten genotypes  
presented as a likelihood ratio. The most likely genotype has a score of 0, and every other genotype's likelihood is a ratio
to this value. 

(Add in example here)

These likelihood ratios can be converted to posterior probabilities using the following equation: 
For a more in depth explanation and proof, see here. (Link)

(Add in equation here)

According to Bayesian statistics, we can either provide this equation with a uniform prior (i.e. our expection is that all outcomes  
are equally likely, or provide a prior ourselves. Using uniform priors would give a probability of observing a heterozygous position  
as 0.6 (6 possible heterozygous genotypes = 0.1 * 6) and of observing a homozygous one only 0.4. We know that observed heterozygosity 
is normally much lower (e.g. x in humans, y in dogs) and using a uniform prior is therefore likely to inflate the number of heterozygous  
calls. To avoid this, we can provide our own estimate of heterozygosity. 

### Step 1: Using ANGSD to get an estimated heterozygosity to use as a prior

One way to obtain aheterozygosity estimate is using the realSFS function in ANGSD to calculate a site frequency spectrum. To do this,  
we need the reference genome for the species (or a closely related one) and optionally an ancestral genome if you want to compute the  
frequency of ancestral homozygous, heterozygous or derived homozygous rather than a folded spectrum (homozygous or heterozygous). 

1. Run angsd to generate a site allele frequency file:
```linux
angsd -i $BAM -r chr10 -anc $REF -ref $REF -C 50 -minQ 20 -minMapQ 30 -dosaf 1 -fold 1 -GL 2 -out ${BAM%\.bam}
```  
   As we are folding the spectrum, we can provide the reference genome for both the -anc and -ref parameters.  

2. Run realSFS to get the sfs estimates:
```linux
realSFS $SAF  > ${BAM%\.bam}.ml
```
3. Compute heterozygosity from the sfs estimates:
The output file from realSFS looks something like the below. To gain an estimate of heterozygosity simply divide the first number  
by the second (double check this) for the folded spectrum.

4. Get bootstrap estimates of the site frequency (optional):
```linux
realSFS  $SAF -bootstrap 1000 > ${BAM%\.bam}_bootstraps.ml
```
### Step 2: Implement the equation in Python to get posterior probabilities

Now we have a heterozygosity estmate to use as a prior in our equation, we can implement it in python with a .glf file as input.  
Below is a breakdown of the script used (full script found here). 

1. Import modules needed (time not essential):
```python
import time
import sys
import csv
import numpy as np
```
2. Define functions to use in the script:  
    a.) Prior calculator
    ```python
    def prior_calculator(heterozygosity):
        "Calculate homozygous and heterozygous priors from input heterozygosity value"
        het = np.divide(heterozygosity, 6.0)
        q = np.subtract(1.0, heterozygosity)
        hom = np.divide(q, 4.0)
        priors = [hom, het, het, het, hom, het, het, hom, het, hom]
        return priors
    ```
    b.) Convert genotype likelihoods from strings to floats:
    ```python
    def float_gls(line):
        "Extract gls and convert from strings to floats"
        gls = [float (i) for i in line[2:12]]
        return gls
    ```
    c.) Return the best genotype (score of 0.0):
    ```python
    def extract_best_genotype(gls):
        "Return genotype with likelihood ratio = 0.0 (best genotype)"
        genotypes = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"] 
        for value in range(0,10):
                if gls[value] == 0.0:
                        best_geno = genotypes[value]
        return best_geno
    ```
    d.) Calculate posterior probabilties for all 10 genotypes or just the one with the highest likelihood:
    ```python
    def gl_calculator(gls, priors = [0.1]*10):      # Default = uniform prior
        "Calculate genotype likelihoods for all genotypes"
        #np.set_printoptions(linewidth=np.inf)
        values = np.power(10, gls)                  # Raise 10 to power of gls
        values = np.multiply(values, priors)        # Multiply by prior
        sum_prob = np.sum(values)
        probabilities = np.divide(values, sum_prob) # Obtain list of probabilities by value/sum_prob
        probabilities = np.round(probabilities, decimals = 5)
        # If just want Pr of best genotype use:
        #best_pr = str(np.amax(probabilities))
        #return best_pr
        # Or if want Pr of all genotypes use:
        list_pr = probabilities.tolist()
        list_pr = [str(i) for i in list_pr]
        return list_pr
    ```
    e.) Append line to .gpf file (genotype probability file):
    ```python
    def append_line(file, line):
        with open(file, "a") as target_file:
                target_file.write(line)
    ```
3.  
    
    
    





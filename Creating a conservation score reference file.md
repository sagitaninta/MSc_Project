# Calculating load scores from pre-computed conservation scores

Short overview of how to take pre-computed phyloP and phastcons scores from UCSC and use these to compute load scores for a species of interest.

## Creating a reference file

### Step 1: Retrieving scores from UCSC

Pre-computed conservation scores are available from the UCSC website for different subsets of species (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/). 
The following will use the placental mammal subset of phastcons and phyloP scores from the 46-way vertebrate alignment:  
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/  
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phastCons46way/placentalMammals/

These folders contain a separate wigFix file for either phyloP or phastcons scores for each position in the alignment, aligned to human chromosomes and scaffolds. 
These can be downloaded to a location of your choice using rsync:
```linux
rsync -avz --progress rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/chrY.phyloP46way.placental.wigFix.gz . # Just scores mapping to chrY of the human genome
rsync -avz --progress rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/*.wigFix.gz . # All score files mapped to the human genome
```

### Step 2: Converting files to bed format

As scores are in a compressed wigFix format, they need converting to a bed format using the wig2bed function from bedops.  
1. Create a list of all the wigFix files you want to convert to bed files and run wig2bed:
```linux
ls -l *wigFix.gz | awk '{ print $NF}' > wig_list.txt
```
2. Use a loop to convert to bed files:  
```linux
while read -r line
do
        gunzip -c $line | wig2bed - > ${line%.wigFix.gz}.bed
done < wig_list.txt
```

### Step 3: Filtering bed files for scores above or below a certain threshold

The next step is to filter the bed files for scores above a threshold (e.g. pre-determined value, top or bottom 5% of values etc).
1. Create a list of all the bed files you want to filter:
```linux
ls -l *.bed | awk '{ print $NF}' > bed_list.txt
```
2. Filter out scores of interest using awk:
```linux
while read -r line 
do
        awk '{ if ( $5 >= 0.43 ) { print } }' $line > $filtered_{line}  # Change threshold score as required (I've used 0.43 for phastcons and 1.5 for phylop)
done < bed_list.txt
```

### Step 4: Lifting over scores to the genome of interest

As scores are for an alignment to the human genome, these need to be lifted over to the reference genome of the species of interest.  
This can be done using the UCSC liftOver tool or using the R package rtracklayer. Chain files for lifting over co-ordinates from hg19  
to other species can be found at http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/

1. Download the chain file for the species of interest (in this case the dog reference genome):
```linux
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToCanFam3.over.chain.gz . 
```
2a. Run command line liftOver:

OR  

2b. Run R implementation of liftOver:
```linux
# Loop to run R script on each bed file:
while read -r line 
do
        Rscript liftOver.R  $line 
done < bed_list.txt
```
liftOver.R content:
```R
# Load rtracklayer and command line args
library(rtracklayer)
args <- commandArgs(TRUE)

# Import chain file:
dog_chain <- import.chain("hg19ToCanFam3.over.chain")

# Import bed file with hg19 co-ordinates to be lifted over:
sample_bed <- import.bed(args[1])

# Liftover from hg19 to new genome co-ordinates:
liftover_out <- liftOver(sample_bed, dog_chain)

# Unlist to go from GRangesList to GRanges object:
test_out_grange <- unlist(liftover_out)

# Set name for out file
out_file <- paste("CanFam_",args[1], sep = "")

# Export the GRange object to bed
export(test_out_grange, out_file, format = "bed")
```
### Step 5: Merging bed files and sorting
Now the bed files with scores have been lifted over. these can be merged and sorted using bedops:
```linux
bedops --everything CanFam_filtered_chr* | sort-bed - --max-mem 10G > all_scores_sorted.bed
```

### Step 6: Filtering out duplicates

The liftover process can lead to multiple co-ordinates and therefore scores mapping to the same position so these need  
removing. 
1. Create a list of duplicate positions (keep only first 3 columns of the bed file):
```linux
awk '{print $1 "\t" $2 "\t" $3}' all_scores_sorted.bed | uniq -d > dup_positions.txt
```
2. Use bedops to remove duplicate positions using the dup_positions.txt file created:
```linux
bedops -n all_scores_sorted.bed dup_positions.txt > filtered_all_scores_sorted.bed
```

### Step 7: Intersecting phyloP scores, phastcons scores and ancestral alleles/reference alleles

To be able to use these scores to calculate a load score, we need to know what the ancestral state at each position is. We can either use the reference genome  
alleles or outgroups to determine these (see x for one way to do this). Once you have these in a file, they can be intersected with the bed file containing  
the scores. If using both phyloP and phastcons scores, these can first be overlapped into a single file.
1. Convert fasta file to bed file for co-ordinates of interest (ignore if already have these in a bed file):
```linux
bedtools getfasta -fi AndeanFox_canFam3.1.fa -bed - -bedOut
```
2. Use bedtools to intersect scores and ancestral/reference alleles:
```linux
bedtools intersect -a anc_phylop_phastcons.bed  -b all_phylop_phastcons_scores.bed -wa -wb -sorted | cut -f 1-4,8-9 > anc_scores_ref.bed
```

Now your reference file should be ready to go! 

### References/software used:


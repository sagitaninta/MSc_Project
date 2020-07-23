# Creating a reference file for computing load scores from pre-computed conservation scores

Short overview of how to take pre-computed phyloP and phastcons scores from UCSC and use these to compute load scores for the species of interest.

## Step 1: Retrieving scores from UCSC

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

## Step 2: Converting files to bed format

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

## Step 3: Filtering bed files for scores above or below a certain threshold

The next step is to filter the bed files for scores above a threshold (e.g. pre-determined value, top or bottom 5% of values etc).
This can easily be done using awk:
```linux
while read -r line 
do
        awk '{ if ( $5 >= 0.43 ) { print } }' $line > $filtered_{line}  # Change threshold 
done < bed_list.txt
```


# Read in list of files, compute coverage and extract depth of coverage (column 7) to write to file with sample names:
while read line
do
        COV=$(samtools coverage -r chr1 $line | awk 'FNR==2{print $7}')
        echo -e `basename ${line%.bam}`"\t"${COV} >> modern_sample_coverage.txt
done < modern_sample_list.txt

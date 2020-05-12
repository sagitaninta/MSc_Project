library(rtracklayer)


# Import chain file:
#dog_chain <- import.chain("hg19ToCanFam2.over.chain")
human_chain <- import.chain("hg19ToHg38.over.chain")

# Import bed file with hg19 co-ordinates to be lifted over:
sample_bed <- import.bed("chrY_sample.bed")

# Liftover from hg19 to new genome co-ordinates:
liftover_out <- liftOver(sample_bed, human_chain)

# Unlist to go from GRangesList to GRanges object:
test_out_grange <- unlist(liftover_out)

# Export the GRange object to bed
export(test_out_grange, "test_out.bed")

read.csv("phylop_tests/bed_files.txt")

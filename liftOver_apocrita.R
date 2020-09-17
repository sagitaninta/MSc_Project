library(rtracklayer) # Load R tracklayer

args <- commandArgs(TRUE) # Take arguments from command line

# Import chain file:
dog_chain <- import.chain("hg19ToCanFam3.over.chain")

# Import bed file with hg19 and CanFan3.1 co-ordinates to be lifted over:
sample_bed <- import.bed(args[1]))

# Liftover from hg19 to new genome co-ordinates:
liftover_out <- liftOver(sample_bed, dog_chain)

# Unlist to go from GRangesList to GRanges object:
test_out_grange <- unlist(liftover_out)

# Set out name for file:
out_file <- paste("hg38_",args[1], sep = "")

# Export the GRange object to bed
export(test_out_grange, out_file, format = "bed")


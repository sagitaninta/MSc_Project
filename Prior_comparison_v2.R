## Set working directory and load libraries:
library(ggplot2)
library(wesanderson)

setwd("~/Masters/Project/11_05_post_pr")

## Read in data:

args = commandArgs(trailingOnly = TRUE) # Need is specifying file in command line

tsv <- read.csv(args[1], sep = "\t", header = FALSE)  # If supplying in files as arguments after code
tsv <- read.csv("AL3185_test_values.gpf", sep = "\t", header = FALSE)  # If supplying file name inside code

tsv <- tsv[1:6] # Get rid of empty column 

colnames(tsv) <- c("chr", "position", "genotype", "probability", "zygosity", "prior") # Specify column names

tsv$prior <- relevel(tsv$prior, ref = "uniform")


tsv2 <- tsv

tsv[7] <- "AL3185"
tsv2[7] <- "CK012"
tsv3 <- rbind(tsv, tsv2)

## Create plot:
# Set colours to use:
cols = wes_palette("Zissou1", n = 4, type = "continuous") # Create palette of specified size (will take first 3 colours in palette of 10 colours)

#Extract sample name:
sample <- sub(pattern = "(.*)_.*", x = args[1], replacement = "\\1")  

#Plot graph:
probability_plot <- ggplot(data = tsv, aes(fill = prior, x = probability)) + 
  theme_minimal() +
  #facet_wrap( ~zygosity, nrow = 2, scales = "free_y") +
  #label_both() +
  facet_grid()
  geom_histogram(binwidth = 0.05, position ="dodge") +   # Histogram with bin width specified
  scale_fill_manual(values = cols) +
  ggtitle(paste("Posterior probabilities for", "sample", sep = " ")) + 
  theme(plot.title = element_text(size = 10)) 
  
  
probability_plot

pdf(file = paste(args[2], ".pdf", sep = ""), width = 6, height = 6)
probability_plot
dev.off()

head(test)

pdf(file = "test.pdf", width = 6, height = 6)
probability_plot
dev.off()

test = "sample name"



probability_plot <- ggplot(data = tsv3, aes(fill = prior, x = probability)) + 
  theme_minimal() +
  #facet_wrap( ~zygosity, nrow = 2, scales = "free_y") +
  #label_both() +
  facet_grid(rows = vars(zygosity), cols = vars(V7), scales = "free_y") +
  geom_histogram(binwidth = 0.05, position ="dodge") +   # Histogram with bin width specified
  scale_fill_manual(values = cols) +
  ggtitle(paste("Posterior probabilities for", "sample", sep = " ")) + 
  theme(plot.title = element_text(size = 10)) 

probability_plot

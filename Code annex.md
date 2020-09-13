# Appendix x: Code

## Creating reference files 
### PhyloP and Phastcons

Test

### SIFT


## Genotype likelihoods and posterior probabilities
Posterior probabilities from .glf file:
```python
import time
import sys
import csv
import numpy as np

start = time.time()

###### Functions ########
def prior_calculator(heterozygosity):
	"Calculate hom and het priors from input heterozygosity value"
	#np.set_printoptions(linewidth=np.inf) #Set numpy print options if needed
	het = np.divide(heterozygosity, 6.0)
	q = np.subtract(1.0, heterozygosity)
	hom = np.divide(q, 4.0)
	priors = [hom, het, het, het, hom, het, het, hom, het, hom]
	return priors

def float_gls(line):
	"Extract gls and convert from strings to floats"
	gls = [float (i) for i in line[2:12]]
	return gls

def extract_best_genotype(gls):
	"Return genotype with likelihood ratio = 0.0 (best genotype)"
	genotypes = ["AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"]   #List genotypes alphabetically
	for value in range(0,10):
		if gls[value] == 0.0:
			best_geno = genotypes[value]     # Store genotype with highest gl
	return best_geno

def gl_calculator(gls, priors = [0.1]*10):               # Default = uniform prior
	"Calculate genotype likelihoods for all genotypes"
	#np.set_printoptions(linewidth=np.inf)
	values = np.power(10, gls)            # Raise 10 to power of gls
	values = np.multiply(values, priors)   # Multiply by prior
	sum_prob = np.sum(values)
	probabilities = np.divide(values, sum_prob) # Obtain list of probabilities by value/sum_prob
	probabilities = np.round(probabilities, decimals = 5)
	# If just want Pr of best genotype use:                #Note: Change here for extracting pr of best genotype or all genotypes 
	#best_pr = str(np.amax(probabilities))
	#return best_pr
	# Or if want Pr of all genotypes use:
	list_pr = probabilities.tolist()
	list_pr = [str(i) for i in list_pr]
	return list_pr

def append_line(file, line):
	with open(file, "a") as target_file:
		target_file.write(line)

###### Code  #######
# Specify input files from command line arguments:
glf_file = sys.argv[1]
output_file = sys.argv[2]

# Calculate heterozygosity/input heterozygosity value:
heterozygosity = 0.0014

# Write posterior probabilities to file:
with open(glf_file, "r") as file:
	glf = csv.reader(file, delimiter = "\t")
	priors = prior_calculator(heterozygosity)   # Comment out if using uniform prior
	for line in glf:
		gls = float_gls(line)
		if np.sum(gls) == 0:
			continue
		else:
			post_probabilities = gl_calculator(gls, priors) # Add in priors if not uniform; if uniform just provide gls (default priors = 0.1)
			bed_pos = str(int(line[1])-1) # Add in 0 based co-ordinate to convert to bed file
			out_line = [line[0], bed_pos, line[1]] + post_probabilities  # Note: Change gl_calculator function output to get probabilties for all genotypes or for just the bestgenotype
			out_line = "\t".join(out_line) # Separate values by tabs
			out_line = out_line + "\n" 
			append_line(output_file, out_line)  # Add line to out file
```

## Overall pipeline

## Plots
Heterozygosity plot:
```R
library(ggplot2)
library(wesanderson)

setwd("~/Masters/Project/Heterozygosity_sfs_estimates/")

############## Code for heterozygosity plot with CIs ##############
##### Get heterozygosity values and CIs from bootstraps #####
## Set file lists:
file_list <- readLines("ml_list.txt")
bootstrap_list <- readLines("bootstrap_list.txt")

## Set vectors to store values in:
all_samples <- c()     #Sample names
all_het <- c()         # Heterozygosity values
all_upper <- c()       # Upper CIs
all_lower <- c()       # Lower CIs

## Loop through files line by line to get names, het values and CIs for each sample:
for(line in 1:length(file_list)) {
  # Scan in files:
  sfs <- scan(file_list[line])
  bootstraps <- read.csv(bootstrap_list[line], header = FALSE, sep = " ")
  # Extract sample name and heterozygosity value:
  sample <- sub(pattern = "(.*)_.*", x = file_list[line], replacement = "\\1")  
  het <- sfs[2]/(sum(sfs))
  # Calculate CIs from bootstrap values:
  all_bootstrap_hets <- c()
  for(line in 1:nrow(bootstraps)) {
  # Calculate heterozygosity for each bootstrap estimate
    bootstrap_het <- (bootstraps[line, 2])/(bootstraps[line, 2]+bootstraps[line, 1])
    # Add to list of all bootstrap het estimates for the sample
    all_bootstrap_hets <- append(all_bootstrap_hets, bootstrap_het) 
  }
  # Take mean het value away from each bootstrap het estimate:
  mean_diffs <- all_bootstrap_hets - het 
  # Extract quantiles from mean difference values:
  quantiles <- quantile(mean_diffs, c(0.025, 0.975))  
  # Het value +/- quantiles to get CI:
  ci <- het - c(quantiles[2], quantiles[1]) 
  # Add values to vector of values for all samples:
  all_samples <- append(all_samples, sample)
  all_het <- append(all_het, het)
  all_upper <- append(all_upper, ci[1])
  all_lower <- append(all_lower, ci[2])
}

###### Create dataframe of values to use for graph #####
#Set group if wanted for key:
group <- c("Dog - Breed", "Dog - Breed", "Wolf", "Dog - Non-breed", "Dog - Non-breed", 
	"Wolf", "Coyote", "Wolf", "Wolf", "Wolf")

# Set names if want different names to sample names in the graph:
sample_names <- c("LabradorRetriever02", "KingCharlesCavalier02", "NorwegianWolf01",
	"ChineseIndigenousDog05", "VietnameseIndigenousDog04", "ChineseWolf05", 
	"CalifornianCoyote", "IranianWolf", "PortugeseWolf", "YellowstoneWolf02")

## Create dataframe of samples, heterozygosity values and CIs:
het_df <- data.frame(all_samples, sample_names, I(group), all_het, all_lower, all_upper)
colnames(het_df) <- c("Sample", "Names", "Group", "Heterozygosity", "LowerCI", "UpperCI")

## Order dataframe by group:
#het_df <- het_df[order(het_df$Group),] # Sorts alphabetically
het_df <- het_df[c(1:2, 4:5, 3, 6, 8:10, 7),] # Sort manually

# Change sample and group to factors (so orders and colours by factor not per bar):
het_df$Sample <- factor(het_df$Sample, levels = het_df$Sample) # If using original names
het_df$Names <- factor(het_df$Names, levels = het_df$Names) # If using shortened names
het_df$Group <- factor(het_df$Group, levels = unique(het_df$Group)) # Keeps order of df rows

## Plot graph with ggplot: 
# Set colours:
cols = wes_palette("IsleofDogs1", n = 12, type = "continuous") # Create palette

# Plot graph vertical bars:
heterozygosity_plot <- ggplot(data = het_df, aes(x = Names, y = Heterozygosity, fill = Group)) +
  geom_bar(stat = "Identity") +
  theme_classic() + xlab(NULL) +
  scale_fill_manual(values = cols) + # Specify colours
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) + # Add error bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Sample names vertical
  ggtitle("Heterozygosity estimates for chromosome 10 across ten canine samples") + # Title
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5, # Title position
  margin = margin(10, 10, 10, 10), vjust = 5)) + 
  theme(axis.text.x = element_text(vjust = 0.5)) + # Move sample labels
  theme(axis.title.y = element_text(vjust = 2))  #+ # Move y axis title
  
heterozygosity_plot

##### Save plot as pdf #####
pdf(file = "Heterozygosity_plot_verticalbars.pdf", width = 12, height = 8)
heterozygosity_plot
dev.off()

write.table(het_df, file = "all_het.csv", sep = "\t")
```


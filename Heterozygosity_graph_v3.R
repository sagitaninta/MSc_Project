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
    bootstrap_het <- (bootstraps[line, 2])/(bootstraps[line, 2]+bootstraps[line, 1]) # Calculate heterozygosity for each bootstrap estimate
    all_bootstrap_hets <- append(all_bootstrap_hets, bootstrap_het) # Add to list of all bootstrap het estimates for the sample
  }
  mean_diffs <- all_bootstrap_hets - het   # Take mean het value away from each bootstrap het estimate
  quantiles <- quantile(mean_diffs, c(0.025, 0.975))  ## Extract quantiles from mean difference values
  ci <- het - c(quantiles[2], quantiles[1])  # Het value +/- quantiles to get CI
  # Add values to vector of values for all samples:
  all_samples <- append(all_samples, sample)
  all_het <- append(all_het, het)
  all_upper <- append(all_upper, ci[1])
  all_lower <- append(all_lower, ci[2])
}

###### Create dataframe of values to use for graph #####
#Set group if wanted for key:
group <- c("Dog - Breed", "Dog - Breed", "Wolf", "Dog - Non-breed", "Dog - Non-breed", "Wolf", "Coyote", "Wolf", "Wolf", "Wolf")

# Set names if want different names to sample names in the graph:
sample_names <- c("LabradorRetriever02", "KingCharlesCavalier02", "NorwegianWolf01", "ChineseIndigenousDog05",
           "VietnameseIndigenousDog04", "ChineseWolf05", "CalifornianCoyote", "IranianWolf",
           "PortugeseWolf", "YellowstoneWolf02")

## Create dataframe of samples, heterozygosity values and CIs:
het_df <- data.frame(all_samples, sample_names, I(group), all_het, all_lower, all_upper)
colnames(het_df) <- c("Sample", "Names", "Group", "Heterozygosity", "LowerCI", "UpperCI")

## Order dataframe by group:
#het_df <- het_df[order(het_df$Group),] # Sorts alphabetically
het_df <- het_df[c(1:2, 4:5, 3, 6, 8:10, 7),] # Sort manually

# Change sample and group to factors (so orders and colours by factor not per bar):
het_df$Sample <- factor(het_df$Sample, levels = het_df$Sample) # If using original sample names
het_df$Names <- factor(het_df$Names, levels = het_df$Names) # If using shortened/ edited names
het_df$Group <- factor(het_df$Group, levels = unique(het_df$Group)) # Keeps order of df rows

## Plot graph with ggplot: 
# Set colours:
cols = wes_palette("IsleofDogs1", n = 12, type = "continuous") # Create palette of specified size (will take first 3 colours in palette of 10 colours)

# Plot graph vertical bars:
heterozygosity_plot <- ggplot(data = het_df, aes(x = Names, y = Heterozygosity, fill = Group)) + geom_bar(stat = "Identity") +
  theme_classic() + xlab(NULL) +
  scale_fill_manual(values = cols) + # Specify colours
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1) + # Add error bars
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Sample names vertical
  ggtitle("Heterozygosity estimates for chromosome 10 across ten canine samples") + # Add title 
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(10, 10, 10, 10), vjust = 5)) + # Modify title size and position
  theme(axis.text.x = element_text(vjust = 0.5)) + # Move sample labels
  theme(axis.title.y = element_text(vjust = 2))  #+ # Move y axis title
  
heterozygosity_plot

##### Save plot as pdf #####
pdf(file = "Heterozygosity_plot_verticalbars.pdf", width = 12, height = 8)
heterozygosity_plot
dev.off()

write.table(het_df, file = "all_het.csv", sep = "\t")

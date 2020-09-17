#### Downsampling tests ####

library(ggplot2)
library(wesanderson)
library(cowplot)

setwd("D:/Debs/downsampling_tests/")

auto_phastcons <- read.csv("all_autocons_phastcons.tsv", sep = "\t", header = FALSE)

auto_phylop <- read.csv("all_autocons_phylop.tsv", sep = "\t", header = FALSE)

auto_sift <- read.csv("all_autosift.tsv", sep = "\t", header = FALSE)

colnames(auto_phastcons) <- c("Coverage", "Sample_ID", "Total_positions", "Hom_anc_tv", "Hom_anc", "Hom_tv", "Het", "Phastcons_score", "Phastcons_per_hom")
colnames(auto_phylop) <- c("Coverage", "Sample_ID", "Total_positions", "Hom_anc_tv", "Hom_anc", "Hom_tv", "Het", "Phylop_score", "Phylop_per_hom")
colnames(auto_sift) <- c("Coverage", "Sample_ID", "Total_positions", "Hom_anc_tv", "Hom_anc", "Hom_tv", "Sift_score", "Sift_per_hom")

auto_phastcons <- subset(auto_phastcons, subset = Sample_ID != "LUPWCHN00013") # Remove old sample
auto_phylop <- subset(auto_phylop, subset = Sample_ID != "LUPWCHN00013")
auto_sift <- subset(auto_sift, subset = Sample_ID != "LUPWCHN00013")

## Phastcons plot

colours = wes_palette("FantasticFox1", n = 7, type = "continuous")

phastcons_plot <- ggplot(auto_phastcons, aes(x = Coverage, y = Total_positions)) +
  scale_colour_manual(values = colours) +
  geom_point(aes(colour = Sample_ID)) +
  ggtitle("Phastcons scores") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 1.0, margin = margin(10, 10, 10, 10), vjust = 5)) + # Modify title size and position
  labs(x = "Downsampled coverage", y = "Total positions", color = "Sample ID") +
  theme_classic()

phastcons_plot

all_plots <- plot_grid(phastcons_plot_cov, phastcons_plot_pos, sift_plot, nrow = 3, labels = "AUTO")

summary(lm(Phastcons_per_hom~Total_positions, data = auto_phastcons))

## PhyloP plot

phylop_plot <- ggplot(auto_phylop, aes(x = Coverage, y = Total_positions)) +
  scale_colour_manual(values = colours) +
  geom_point(aes(colour = Sample_ID)) +
  ggtitle("Phylop scores") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 1.0, margin = margin(10, 10, 10, 10), vjust = 5)) + # Modify title size and position
  labs(x = "Downsampled coverage", y = "Total positions", color = "Sample ID") +
  theme_classic()

phylop_plot

## Sift plot

sift_plot <- ggplot(auto_sift, aes(x = Coverage, y = Total_positions)) +
  scale_colour_manual(values = colours) +
  geom_point(aes(colour = Sample_ID)) +
  ggtitle("Sift scores") +
  theme(plot.title = element_text(size= 10, face = "bold", hjust = 1.0, margin = margin(10, 10, 10, 10), vjust = 5)) + # Modify title size and position
  labs(x = "Downsampled coverage", y = "Total positions", color = "Sample ID") +
  theme_classic()

sift_plot

## Combine plots

title <- ggdraw() +
  draw_label("Homozygous transversions vs coverage", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 


all_plots <- plot_grid(phastcons_plot, phylop_plot, sift_plot, nrow = 3, labels = "AUTO")
plot_title <- plot_grid(title, all_plots, ncol = 1, rel_heights = c(0.1,1))
plot_title

pdf(file = "Downsampling_tests - Score per hom, cov.pdf", width = 8, height = 14)
plot_title
dev.off()

##### Analysis script #####

## Set working directory and load packages:

setwd("D:/Debs/Results/")

library(ggplot2)
library(wesanderson)
library(cowplot)
library(ggrepel)
library(ggforce)
library(dplyr)

options(scipen = 999) # Prevent scientific notation on axis

### Read in data and filter out samples to exclude:

## Metadata:
metadata <- read.csv("Metadata - Metadata for samples new v4.tsv", sep = "\t", header = TRUE)

## Phastcons

all_phastcons <- read.csv("all_phastcons_v5.tsv", sep = "\t", header = FALSE)
colnames(all_phastcons) <- c("Sample_ID", "Total_positions", "Hom_anc_tv", "Hom_anc", "Hom_tv", "Het", "Phastcons_score", "Phastcons_per_hom")
phastcons_merge <- merge(metadata, all_phastcons, by = "Sample_ID", incomparables = NA)
phastcons_merge <- subset(phastcons_merge, subset = Sample_ID != "CGG6" & Sample_ID != "LUPWCHN00013" & Sample_ID != "RedWolf_RW3")

## PhyloP
all_phylop <- read.csv("all_phylop_v5.tsv", sep = "\t", header = FALSE)
colnames(all_phylop) <- c("Sample_ID", "Total_positions", "Hom_anc_tv", "Hom_anc", "Hom_tv", "Het", "Phylop_score", "Phylop_per_hom")
phylop_merge <- merge(metadata, all_phylop, by = "Sample_ID", incomparables = NA)
phylop_merge <- subset(phylop_merge, subset = Sample_ID != "CGG6" & Sample_ID != "LUPWCHN00013" & Sample_ID != "RedWolf_RW3")

## SIFT
all_sift <- read.csv("all_sift_v5.tsv", sep = "\t", header = FALSE)
colnames(all_sift) <- c("Sample_ID", "Total_positions", "Hom_anc_tv", "Hom_anc", "Hom_tv", "Sift_score", "Sift_per_hom" & Sample_ID != "rw_texas2")
sift_merge <- merge(metadata, all_sift, by = "Sample_ID", incomparables = NA)
sift_merge <- subset(sift_merge, subset = Sample_ID != "CGG6" & Sample_ID != "LUPWCHN00013" & Sample_ID != "RedWolf_RW3")

############################### Plots - all data #######################################

# Set colours
anc_mod_cols <- c("#F98400", "#B40F20", "#5BBCD6")
#anc_mod_cols <- wes_palette("FantasticFox1", n = 3, type = "continuous")

#### Total positions by coverage #######################

pc_cov_pos <- ggplot(phastcons_merge, aes(x = Coverage, y = Total_positions)) +
  geom_point(aes(colour = Type)) +
  scale_colour_manual(values = anc_mod_cols) +
  geom_vline(xintercept = 4, linetype = "dashed") +
  #geom_hline(yintercept = 105000000, linetype = "dotted") +
  labs(x = "Coverage", y = "Total positions") +
  ggtitle("Phastcons scores") +
  theme_classic()
pc_cov_pos

py_cov_pos <- ggplot(phylop_merge, aes(x = Coverage, y = Total_positions)) +
  geom_point(aes(colour = Type)) +
  scale_colour_manual(values = anc_mod_cols) +
  geom_vline(xintercept = 4, linetype = "dashed") +
  #geom_hline(yintercept = 85000000, linetype = "dotted") +
  labs(x = "Coverage", y = "Total positions") +
  ggtitle("Phylop scores") +
  theme_classic()
py_cov_pos

si_cov_pos <- ggplot(sift_merge, aes(x = Coverage, y = Total_positions)) +
  geom_point(aes(colour = Type)) +
  scale_colour_manual(values = anc_mod_cols) +
  geom_vline(xintercept = 4, linetype = "dashed") +
  #geom_hline(yintercept = 14500000, linetype = "dotted") +
  labs(x = "Coverage", y = "Total positions") +
  ggtitle("Sift scores") +
  theme_classic()
si_cov_pos

# Combine plots:
title <- ggdraw() +
  draw_label("Total positions by coverage", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

plots_cov_pos <- plot_grid(pc_cov_pos, py_cov_pos, si_cov_pos, nrow = 3, labels = "AUTO")
plots_cov_pos <- plot_grid(title, plots_cov_pos, ncol = 1, rel_heights = c(0.1,1))
plots_cov_pos

pdf(file = "Thesis figs - Total positions by cov v1.pdf", width = 8, height = 14)
plots_cov_pos
dev.off()

#### Score by coverage ##############################

# Phastcons:
pc_cov_score <- ggplot(phastcons_merge, aes(x = Coverage, y = Phastcons_per_hom)) +
  scale_colour_manual(values = anc_mod_cols) +
  geom_point(aes(colour = Type)) +
  #scale_colour_manual(values = anc_mod_cols) +
  geom_vline(xintercept = 4, linetype = "dashed") +
  labs(x = "Coverage", y = "Score per homozygous position") +
  ggtitle("Phastcons scores") +
  theme_classic()
pc_cov_score

# Phylop:
py_cov_score <- ggplot(phylop_merge, aes(x = Coverage, y = Phylop_per_hom)) +
  scale_colour_manual(values = anc_mod_cols) +
  geom_point(aes(colour = Type)) +
  geom_vline(xintercept = 4, linetype = "dashed") +
  labs(x = "Coverage", y = "Score per homozygous position") +
  ggtitle("Phylop scores") +
  theme_classic()
py_cov_score

# Sift:
si_cov_score <- ggplot(sift_merge, aes(x = Coverage, y = Sift_per_hom)) +
  geom_point(aes(colour = Type)) +
  scale_colour_manual(values = anc_mod_cols) +
  geom_vline(xintercept = 4, linetype = "dashed") +
  labs(x = "Coverage", y = "Score per homozygous position") +
  ggtitle("Sift scores") +
  theme_classic()
si_cov_score

# Combine plots:
title <- ggdraw() +
  draw_label("Scores by coverage", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

plots_cov_score <- plot_grid(pc_cov_score, py_cov_score, si_cov_score, nrow = 3, labels = "AUTO")
plot_cov_score <- plot_grid(title, plots_cov_score, ncol = 1, rel_heights = c(0.1,1))
plot_cov_score

pdf(file = "Thesis figs - scores by coverage.pdf", width = 8, height = 14)
plot_cov_score
dev.off()


###  Stats

glm1 <- lm(Phastcons_per_hom~Coverage+Epoch, data = phastcons_subset)
summary(glm1)


################################# Plots - Filtered data #############################################
#### Subset data:
# Phastcons:
phastcons_subset <- subset(phastcons_merge, subset = Coverage >= 4)

# PhyloP:
phylop_subset <- subset(phylop_merge, subset = Coverage >= 4)

# Sift:
sift_subset <- subset(sift_merge, subset = Coverage >= 4)
### Set colours:

#group_cols = wes_palette("FantasticFox1", n = 4)

# Merge datasets

all_merge <- merge(phastcons_subset, phylop_subset, by = "Sample_ID", incomparables = NA)
all_merge <- merge(all_merge, sift_subset, by = "Sample_ID", incomparables = NA)
all_merge_meta <- all_merge[1:13] 
  
Phastcons_per_hom <- all_merge$Phastcons_per_hom
Phylop_per_hom <- all_merge$Phylop_per_hom
Sift_per_hom <- all_merge$Sift_per_hom
meta <- all_merge[,1:13]

all_data <- cbind(meta, Phastcons_per_hom, Phylop_per_hom, Sift_per_hom)

#### Scores by group ###############################

anc_mod_bw <- c("chocolate", "darkgrey", "darkslategrey")
#group_cols = wes_palette()

# Phastcons:
pc_group <- ggplot(phastcons_subset, aes(x = Epoch, y = Phastcons_per_hom, label = Graph_name)) +
  geom_violin(alpha = 0.4, colour = "darkgrey", aes(fill = Epoch), draw_quantiles = c(0.5)) +
  geom_jitter(aes(colour = Type), width = 0.2, shape = 4, size = 3) +
  geom_text_repel(data = phastcons_subset %>% filter((Phastcons_per_hom > 0.000175 | Graph_name == "villdog_vietnam" |
                  Graph_name == "villdog_china" | Graph_name == "villdog_portugal" | 
                    Graph_name == "din_australia")), size = 3) +
  scale_colour_manual(values = anc_mod_bw) +
  ggtitle("Phastcons score") +
  labs(x = "Group", y = "Phastcons score per homozygous position") +
  scale_x_discrete(limits=c("Holocene_America", "Holocene_Eur_Mid_East", "Holocene_East",  "Pleistocene", "Holocene_dog", "Holocene_other"),
                   labels = c("America","Europe & Middle East", "East Asia",  "Pleistocene", "Dogs", "Other canids")) +
  guides(fill = FALSE) +
  #legend() +
  theme_classic()
pc_group

# Phylop:
py_group <- ggplot(phylop_subset, aes(x = Epoch, y = Phylop_per_hom, label = Graph_name)) +
  geom_violin(alpha = 0.4, colour = "darkgrey", aes(fill = Epoch), draw_quantiles = c(0.5)) +
  geom_jitter(aes(colour = Type), width = 0.2, shape = 4, size = 3) +
  geom_text_repel(data = phylop_subset %>% filter((Phylop_per_hom > 0.00027 | Graph_name == "villdog_vietnam" |
                                                     Graph_name == "villdog_china" | Graph_name == "villdog_portugal" |
                                                     Graph_name == "din_australia")), size = 3) +
  scale_colour_manual(values = anc_mod_bw) +
  ggtitle("Phylop score") +
  labs(x = "Group", y = "Phylop score per homozygous position") +
  scale_x_discrete(limits=c("Holocene_America", "Holocene_Eur_Mid_East", "Holocene_East",  "Pleistocene", "Holocene_dog", "Holocene_other"),
                   labels = c("America","Europe & Middle East", "East Asia",  "Pleistocene", "Dogs", "Other canids")) +
  guides(fill = FALSE) +
  theme_classic()
py_group

# Sift:
si_group <- ggplot(sift_subset, aes(x = Epoch, y = Sift_per_hom, label = Graph_name)) +
  geom_violin(alpha = 0.4, colour = "darkgrey", aes(fill = Epoch), draw_quantiles = c(0.5)) +
  geom_jitter(aes(colour = Type), width = 0.2, shape = 4, size = 3) +
  geom_text_repel(data = sift_subset %>% filter((Sift_per_hom > 0.00004|  Graph_name == "villdog_vietnam" |
                                                   Graph_name == "villdog_china" | Graph_name == "villdog_portugal" |
                                                   Graph_name == "din_australia" | Graph_name == "rw_alabama_his")), size = 3) +
  scale_colour_manual(values = anc_mod_bw) +
  ggtitle("Sift score") +
  labs(x = "Group", y = "Sift score per homozygous position") +
  scale_x_discrete(limits=c("Holocene_America", "Holocene_Eur_Mid_East", "Holocene_East",  "Pleistocene", "Holocene_dog", "Holocene_other"),
                   labels = c("America","Europe & Middle East", "East Asia",  "Pleistocene", "Dogs", "Other canids")) +
  guides(fill = FALSE) +
  theme_classic()
si_group

# Combine plots:
title <- ggdraw() +
  draw_label("Load scores by region", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

plots_date_group <- plot_grid(pc_group, py_group, si_group, nrow = 3, labels = "AUTO")
plots_date_group <- plot_grid(title, plots_date_group, ncol = 1, rel_heights = c(0.1,1))
plots_date_group

pdf(file = "Thesis figs - plots by region main.pdf", width = 8, height = 14)
plots_date_group
dev.off()

## Stats

lm1 <- lm(Phastcons_per_hom~Epoch+Coverage, data = phastcons_subset)
anova(lm1)
summary(lm1)

#### Scores by score ##########################
group_cols = wes_palette("FantasticFox1", n = 7, type = "continuous")

# Phastcons/Phylop:
pc_py <- ggplot(phastcons_subset, aes(x = Phastcons_per_hom, y = phylop_subset$Phylop_per_hom)) +
  geom_point(aes(colour = Epoch, shape = Type)) +
  scale_colour_manual(values = group_cols, labels = c("Holocene - America", "Holocene - Dogs", "Holocene - Eurasia", "Holocene - Other canids", "Pleistocene - Siberia")) +
  ggtitle("Phylop / Phastcons") +
  labs(x = "Phastcons per homozygous position", y = "Phylop per homozygous position", colour = "Group", shape = "Sample type") +
  theme_classic()
pc_py

# Phastcons/Sift:
pc_si <- ggplot(phastcons_subset, aes(x = Phastcons_per_hom, y = sift_subset$Sift_per_hom)) +
  geom_point(aes(colour = Epoch, shape = Type), ) +
  scale_colour_manual(values = group_cols, labels = c("Holocene - America", "Holocene - Dogs", "Holocene - Eurasia", "Holocene - Other canids", "Pleistocene - Siberia")) +
  ggtitle("Sift / Phastcons") +
  labs(x = "Phastcons per homozygous position", y = "Sift per homozygous position", colour = "Group", shape = "Sample type") +
  theme_classic()
pc_si

# Phylop/Sift:
py_si <- ggplot(phylop_subset, aes(x = Phylop_per_hom, y = sift_subset$Sift_per_hom)) +
  geom_point(aes(colour = Epoch, shape = Type)) +
  scale_colour_manual(values = group_cols, labels = c("Holocene - America", "Holocene - Dogs", "Holocene - Eurasia", "Holocene - Other canids", "Pleistocene - Siberia")) +
  ggtitle("Phylop / Sift") +
  labs(x = "Phylop per homozygous position", y = "Sift per homozygous position") +
  theme_classic()
py_si

# Combine plots:
title <- ggdraw() +
  draw_label("Correlation between load scores", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

plots_score_correl <- plot_grid(pc_py, pc_si, py_si, nrow = 3, labels = "AUTO")
plots_score_correl <- plot_grid(title, plots_date_age, ncol = 1, rel_heights = c(0.1,1))
plots_score_correl

pdf(file = "Thesis figs - Correl between scores.pdf", width = 8, height = 14)
plots_score_correl
dev.off()

## Stats

cor.test(phastcons_subset$Phastcons_per_hom, phylop_subset$Phylop_per_hom, alternative = "greater")
cor.test(phastcons_subset$Phastcons_per_hom, sift_subset$Sift_per_hom, alternative = "greater")
cor.test(phylop_subset$Phylop_per_hom, sift_subset$Sift_per_hom, alternative = "greater")

################################# Plots - Regions ##############################################

#### Subset data:
# Phastcons:
phastcons_subset_america <- subset(phastcons_subset, subset = Group == "North_America")
phastcons_subset_europe <- subset(phastcons_subset, subset = Group == "Europe" | Region_grouping == "MiddleEast_wolf")
phastcons_subset_asia <- subset(phastcons_subset, subset = Group == "Asia" & Region_grouping != "MiddleEast_wolf")

# Phylop:
phylop_subset_america <- subset(phylop_subset, subset = Group == "North_America")
phylop_subset_europe <- subset(phylop_subset, subset = Group == "Europe" | Region_grouping == "MiddleEast_wolf")
phylop_subset_asia <- subset(phylop_subset, subset = Group == "Asia" & Region_grouping != "MiddleEast_wolf")

# Sift:
sift_subset_america <- subset(sift_subset, subset = Group == "North_America")
sift_subset_europe <- subset(sift_subset, subset = Group == "Europe" | Region_grouping == "MiddleEast_wolf")
sift_subset_asia <- subset(sift_subset, subset = Group == "Asia" & Region_grouping != "MiddleEast_wolf")

### Set colours:
am_cols <- wes_palette("Darjeeling1", n = 10, type = "continuous")

#### North America #################################

# Phastcons
pc_region_am <- ggplot(phastcons_subset_america, aes(x = Region_grouping, y = Phastcons_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping), na.rm = FALSE) +
  geom_text_repel(data = phastcons_subset_america %>% filter(Graph_name == "stlawrenceisland" | Graph_name == "saskatchewan" |
                  Graph_name == "greenland" | Graph_name == "isleroyale" | Graph_name =="montana" | Graph_name == "mexico" |  Region_grouping == "Red_wolf" | Graph_name == "coy_wyoming" |
                    Graph_name == "coy_alabama" | Region_grouping == "Hybrid"), size = 3, point.padding = 0.1) +
  #geom_text_repel(size =3) +
  scale_colour_manual(values = am_cols, name = "Group", limits=c("Alaska_wolf", "Artic_wolf", "Canada_wolf", "Other_wolf", "Hybrid", "Red_wolf", "Coyote"), 
                      labels = c("Alaskan wolves", "Arctic wolves", "Canadian wolves", "US and Mexican wolves", "Hybrids", "Red wolves", "Coyotes")) +
  ggtitle("Phastcons score") +
  labs(x = "Group", y = "Score per homozygous position") +
  scale_x_discrete(limits=c("Alaska_wolf", "Artic_wolf", "Canada_wolf", "Other_wolf", "Hybrid", "Red_wolf", "Coyote"),
                   labels = c("Alaska", "Arctic", "Canada", "Other", "Hybrids", "Red wolves", "Coyotes")) +
  guides(fill = FALSE) 
  
pc_region_am


# PhyloP
py_region_am <- ggplot(phylop_subset_america, aes(x = Region_grouping, y = Phylop_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping), na.rm = FALSE) +
  geom_text_repel(data = phylop_subset_america %>% filter(Graph_name == "stlawrenceisland" | Graph_name == "saskatchewan" |
                                                              Graph_name == "greenland" | Graph_name == "isleroyale" | Graph_name =="montana" | Graph_name == "mexico" | Region_grouping == "Red_wolf" | 
                                                              Graph_name == "coy_alabama" | Graph_name == "coy_wyoming" | Region_grouping == "Hybrid"), size = 3, point.padding = 0.1) +
  #geom_text_repel(size =3) +
  scale_colour_manual(values = am_cols, name = "Group", limits=c("Alaska_wolf", "Artic_wolf", "Canada_wolf", "Other_wolf", "Hybrid", "Red_wolf", "Coyote"), 
                      labels = c("Alaskan wolves", "Arctic wolves", "Canadian wolves", "US and Mexican wolves", "Hybrids", "Red wolves", "Coyotes")) +
  ggtitle("Phylop score") +
  labs(x = "Group", y = "Score per homozygous position") +
  scale_x_discrete(limits=c("Alaska_wolf", "Artic_wolf", "Canada_wolf", "Other_wolf", "Hybrid", "Red_wolf", "Coyote"),
                   labels = c("Alaska", "Arctic", "Canada", "Other", "Hybrids", "Red wolves", "Coyotes")) +
  guides(fill = FALSE) 
  
py_region_am

# Sift
si_region_am <- ggplot(sift_subset_america, aes(x = Region_grouping, y = Sift_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping), na.rm = FALSE) +
  geom_text_repel(data = sift_subset_america %>% filter(Graph_name == "isleroyale" | Graph_name =="montana" | Graph_name == "mexico" | Graph_name == "qamanirjuaq" |
                  Graph_name == "greenland" | Region_grouping == "Red_wolf" | Graph_name == "coy_alabama" | Graph_name == "coy_wyoming" |
                    Region_grouping == "Hybrid" | Graph_name == "hy_greatlakes" | Graph_name == "hy_greatplains_his"), size = 3, point.padding = 0.1) +
  #geom_text_repel(size =3) +
  scale_colour_manual(values = am_cols, name = "Group", limits=c("Alaska_wolf", "Artic_wolf", "Canada_wolf", "Other_wolf", "Hybrid", "Red_wolf", "Coyote"), 
                      labels = c("Alaskan wolves", "Arctic wolves", "Canadian wolves", "US and Mexican wolves", "Hybrids", "Red wolves", "Coyotes")) +
  ggtitle("Sift score") +
  labs(x = "Group", y = "Score per homozygous position") +
  scale_x_discrete(limits=c("Alaska_wolf", "Artic_wolf", "Canada_wolf", "Other_wolf", "Hybrid", "Red_wolf", "Coyote"),
                   labels = c("Alaska", "Arctic", "Canada", "Other", "Hybrids", "Red wolves", "Coyotes")) +
  guides(fill = FALSE) 
  
si_region_am

## Combine plots:

title <- ggdraw() +
  draw_label("Scores by region - North America", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

all_plots_am <- plot_grid(pc_region_am, py_region_am, si_region_am, nrow = 3, labels = "AUTO")
all_plots_am <- plot_grid(title, all_plots_am, ncol = 1, rel_heights = c(0.1,1))
all_plots_am

pdf("Thesis figs = All NAm samples main.pdf", width = 8, height = 14)
all_plots_am
dev.off()

#### Europe + Middle East #######################################

## Set colours:
eur_cols <- wes_palette("Zissou1", n = 6, type = "continuous")

# Phastcons
pc_region_eur <- ggplot(phastcons_subset_europe, aes(x = Region_grouping, y = Phastcons_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping), na.rm = FALSE) +
  #geom_text_repel(data = phastcons_subset_europe %>% filter(Region_grouping == "Southeur_wolf" | Region_grouping == "MiddleEast_wolf" |
  #               Region_grouping == "Ancient_wolf" | Graph_name == "dog_ireland4.8k" | Graph_name == "villdog_portugal"), size = 3, point.padding = 0.1) +
  geom_text_repel(size =3) +
  scale_colour_manual(values = eur_cols, name = "Group", limits=c("Southeur_wolf", "Greek_wolf", "Scandinavian_wolf", "MiddleEast_wolf", "Ancient_wolf", "Dog"), 
                      labels = c("Southern Europe wolves", "Greek wolves", "Scandinavian wolves", "Middle Eastern wolves", "Ancient wolves", "Dogs")) +
  ggtitle("Phastcons score") +
  labs(x = "Group", y = "Score per homozygous position") +
  scale_x_discrete(limits=c("Southeur_wolf", "Greek_wolf", "Scandinavian_wolf", "MiddleEast_wolf", "Ancient_wolf", "Dog"),
                   labels = c("South Europe", "Greece", "Scandinavia", "Middle East", "Ancient", "Dogs")) +
  guides(fill = FALSE) 

pc_region_eur

# PhyloP
py_region_eur <- ggplot(phylop_subset_europe, aes(x = Region_grouping, y = Phylop_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping), na.rm = FALSE) +
  #geom_text_repel(data = phylop_subset_europe %>% filter(Region_grouping == "Southeur_wolf" | Region_grouping == "MiddleEast_wolf" |
  #                                                            Region_grouping == "Ancient_wolf" | Graph_name == "dog_ireland4.8k" | Graph_name == "villdog_portugal"), size = 3, point.padding = 0.1) +
  geom_text_repel(size =3) +
  scale_colour_manual(values = eur_cols, name = "Group", limits=c("Southeur_wolf", "Greek_wolf", "Scandinavian_wolf", "MiddleEast_wolf", "Ancient_wolf", "Dog"), 
                      labels = c("Southern Europe wolves", "Greek wolves", "Scandinavian wolves", "Middle Eastern wolves", "Ancient wolves", "Dogs")) +
  ggtitle("Phylop score") +
  labs(x = "Group", y = "Score per homozygous position") +
  scale_x_discrete(limits=c("Southeur_wolf", "Greek_wolf", "Scandinavian_wolf", "MiddleEast_wolf", "Ancient_wolf", "Dog"),
                   labels = c("South Europe", "Greece", "Scandinavia", "Middle East", "Ancient", "Dogs")) +
  guides(fill = FALSE) 

py_region_eur

# Sift
si_region_eur <- ggplot(sift_subset_europe, aes(x = Region_grouping, y = Sift_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping), na.rm = FALSE) +
  #geom_text_repel(data = sift_subset_europe %>% filter(Region_grouping == "Southeur_wolf" | Graph_name == "israel" |
  #                                                            Region_grouping == "Ancient_wolf" | Graph_name == "dog_ireland4.8k" | Graph_name == "villdog_portugal"), size = 3, point.padding = 0.1) +
  geom_text_repel(size =3) +
  scale_colour_manual(values = eur_cols, name = "Group", limits=c("Southeur_wolf", "Greek_wolf", "Scandinavian_wolf", "MiddleEast_wolf", "Ancient_wolf", "Dog"), 
                      labels = c("Southern Europe wolves", "Greek wolves", "Scandinavian wolves", "Middle Eastern wolves", "Ancient wolves", "Dogs")) +
  ggtitle("Sift score") +
  labs(x = "Group", y = "Score per homozygous position") +
  scale_x_discrete(limits=c("Southeur_wolf", "Greek_wolf", "Scandinavian_wolf", "MiddleEast_wolf", "Ancient_wolf", "Dog"),
                   labels = c("South Europe", "Greece", "Scandinavia", "Middle East", "Ancient", "Dogs")) +
  guides(fill = FALSE) 

si_region_eur
## Combine plots:

title <- ggdraw() +
  draw_label("Scores by region - Europe", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

all_plots_eur <- plot_grid(pc_region_eur, py_region_eur, si_region_eur, nrow = 3, labels = "AUTO")
all_plots_eur <- plot_grid(title, all_plots_eur, ncol = 1, rel_heights = c(0.1,1))

pdf("Thesis figs - All Eur Mid East all labels annex.pdf", width = 8, height = 16)
all_plots_eur
dev.off()

### East Asia ##########################################
## Set colours
asia_cols <- wes_palette("Zissou1", n = 7, type = "continuous")

# Phastcons
pc_region_asia <- ggplot(phastcons_subset_asia, aes(x = Region_grouping, y = Phastcons_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping)) +
  scale_colour_manual(values = asia_cols) +
  ggtitle("Phastcons score") +
  geom_text_repel(size = 3) +
  scale_x_discrete(limits=c("Mongolia_wolf", "Tibet_wolf", "Xinjiang_wolf", "China_wolf", "Other_wolf", "Dog"),
                   labels = c("Mongolia", "Tibet", "Xinjiang", "Other China", "India", "Dogs")) +
  guides(fill = FALSE) 
pc_region_asia


# PhyloP
py_region_asia <- ggplot(phylop_subset_asia, aes(x = Region_grouping, y = Phylop_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping)) +
  scale_colour_manual(values = asia_cols) +
  geom_text_repel(size = 3) +
  ggtitle("Phylop score") +
  scale_x_discrete(limits=c("Mongolia_wolf", "Tibet_wolf", "Xinjiang_wolf", "China_wolf", "Other_wolf", "Dog"),
                   labels = c("Mongolia", "Tibet", "Xinjiang", "Other China", "India", "Dogs")) +
  guides(fill = FALSE) 
py_region_asia

# Sift
si_region_asia <- ggplot(sift_subset_asia, aes(x = Region_grouping, y = Sift_per_hom, label = Graph_name)) +
  theme_minimal_grid() +
  theme(axis.text = element_text(size = 10), plot.title = element_text(size = 12), 
        axis.title = element_text(size = 12), legend.title = element_text(size = 12), legend.text = element_text(size = 10)) +
  geom_point(aes(colour = Region_grouping)) +
  scale_colour_manual(values = asia_cols) +
  ggtitle("Sift score", aes(size = 12)) +
  geom_text_repel(size = 3) +
  scale_x_discrete(limits=c("Mongolia_wolf", "Tibet_wolf", "Xinjiang_wolf", "China_wolf", "Other_wolf", "Dog"),
                   labels = c("Mongolia", "Tibet", "Xinjiang", "Other China", "India", "Dogs")) +
  guides(fill = FALSE) 
si_region_asia

## Combine plots:

title <- ggdraw() +
  draw_label("Scores by region - Asia", fontface = "bold", hjust = 0.5) +
  theme(plot.margin = margin(0,0,0,7)) 

all_plots_asia <- plot_grid(pc_region_asia, py_region_asia, si_region_asia, nrow = 3, labels = "AUTO")
all_plots_asia <- plot_grid(title, all_plots_asia, ncol = 1, rel_heights = c(0.1,1))

pdf("Thesis figures - Region Asia all labels.pdf", width = 8, height = 14)
all_plots_asia
dev.off()

#### Map of samples ##################################################################################
# Read in packages
library(plotly)
library(readr)

Sys.setenv("plotly_username"="debgreer") #set up plotly account here https://plot.ly/accounts/login/?action=login 
Sys.setenv("plotly_api_key"="cZTJahK6KAeRJhwmz4Kf")

#Extract metadata for filtered wolf samples (Grey wolf and red wolf:
wolves <- subset(metadata, subset = Epoch == "Holocene_East" | Epoch == "Holocene_Eur_Mid_East" | Epoch == "Holocene_America" | Epoch == "Pleistocene")
wolves <- subset(wolves, subset = Coverage > 4)
wolves <- subset(wolves, subset = Sample_ID != "CGG6" & Sample_ID != "LUPWCHN00013" & Sample_ID != "RedWolf_RW3")

# Set colours:
plotcols <- c("#F21A00", "#E1AF00", "#EBCC2A",  "#78B7C5", "#3B9AB2")

# Plot graph
graph <- list(showframe = FALSE, # Draws a frame around the map
          coastlinecolor = toRGB("white"),
          showland = TRUE,
          landcolor = toRGB("gray80"),
          showcountries = FALSE,
          countrycolor = toRGB("white"), # Lines between countries
          countrywidth = 0.2,
          projection = list(type = 'Mercator')) # Shape of the map go to projection type in https://plotly.com/r/reference/  for list of options
#This section plots onto the map 
wolfplot <- plot_geo(wolves) %>% #Samples refers to above where you read in the csv file 
  add_markers(x = ~Longitude, 
              y = ~Latitude, 
              #size = ~Individuals, # Can change the size of the marker
              symbol = ~Type, # Split by symbol
              symbols = c("square", "circle", "triangle-up"),
              color = ~Analysis_date, # Scale by colour 
              #size = 10, # Should be in marker bracket
              colors = plotcols, # List of colours you want to use for your markers
              opacity=0.70,
              marker = list(
                line = list(color = "black", width = 0.5), 
                size = 12), # 
              #size = -Samples$`Number in country`), change size to no. samples
              text = ~paste(Sample_ID, Graph_name), #put in here what you want to see on hover 
              #hoverinfo = "text" #this then calls what you specified above 
  ) %>%
  #this section draws the underlying map 
  layout(geo = graph, #calls g from above 
         title = '<b> Distribution of Grey wolf and Red wolf samples </b>',
         colorscale = list(title=list(text="Analysis date"))
         ,legend=list(title=list(text='Sample type')) #I have this turned off but included it incase anyone needs it 
  ) 

wolfplot

htmlwidgets::saveWidget(fig, file="Wolf map v1.html") #html file to local computer
api_create(fig, filename = "map_test2") #opens on plotly online 

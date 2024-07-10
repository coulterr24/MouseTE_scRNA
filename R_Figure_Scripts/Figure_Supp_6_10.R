
#### Figure Supp 6 and 10: Stem and Cancer Markers ####

setwd("/Users/Coulter/Desktop/PhD Research/Seurat/Mouse Samples/Mouse_TE_Cancer")

#### Packages Load ####

library(dplyr)
library(patchwork)
library(Seurat)
library(harmony)
library(ggplot2)
library(cowplot)
library(SoupX)
library(DoubletFinder)
library(data.table)
library(parallel)
library(tidyverse)
library(SoupX)
library(ggrepel)

library(ggplot2)
library(gplots)
library(RColorBrewer)
library(viridisLite)
library(Polychrome)
library(circlize)
library(NatParksPalettes)

library(monocle3)
library(ComplexHeatmap)
library(ggExtra)
library(gridExtra)
library(egg)

library(scales)

#### Distal Epithelial and Pseudotime Dataset ####


Epi_Filter <- readRDS(file = "../Distal/20220817_Distal_Epi_Cells.rds" , refhook =  NULL)

Epi_Named <- RenameIdents(Epi_Filter, 
                          '0' = "Spdef+ Secretory", 
                          '1' = "Slc1a3+ Stem/Progenitor", 
                          '2' = "Cebpdhigh/Foxj1- Progenitor",
                          '3' = "Ciliated 1", 
                          '4' = "Ciliated 2", 
                          '5' = "Pax8low/Prom1+ Cilia-forming", 
                          '6' = "Fibroblast-like",
                          '7' = "Slc1a3med/Sox9+ Cilia-forming",
                          '8' = "Selenop+/Gstm2high Secretory")

Epi_Named@active.ident <- factor(x = Epi_Named@active.ident, levels = c( c("Slc1a3+ Stem/Progenitor",
                                                                           "Cebpdhigh/Foxj1- Progenitor",
                                                                           "Slc1a3med/Sox9+ Cilia-forming",
                                                                           "Pax8low/Prom1+ Cilia-forming", 
                                                                           "Fibroblast-like",
                                                                           "Spdef+ Secretory",
                                                                           "Selenop+/Gstm2high Secretory",
                                                                           "Ciliated 1",
                                                                           "Ciliated 2")))

cds <- readRDS(file = "20221101_Distal_Epi_PHATE_Monocle3.rds" , refhook = NULL)

#### Figure Supp 4: Stem Dot Plot ####


stem_features <- c("Krt5","Krt17","Cd44","Prom1","Kit","Aldh1a1","Aldh1a2","Aldh1a3",
                   "Efnb1","Ephb1","Trp63","Sox2","Sox9","Klf4","Rnf43","Foxm1",
                   "Pax8","Nanog","Itga6","Psca","Tcf3","Tcf4","Nrp1","Slc1a3","Tnfrsf19",
                   "Smo","Lrig1","Ezh2","Egr1","Tacstd2","Dusp1","Slc38a2","Malat1",
                   "Btg2","Cdkn1c","Pdk4","Nedd9","Fos","Jun","Junb","Zfp36",
                   "Neat1","Gadd45g","Gadd45b")


stem_dp <- DotPlot(object = Epi_Named,                    # Seurat object
                   assay = 'RNA',                        # Name of assay to use.  Default is the active assay
                   features = stem_features,                  # List of features (select one from above or create a new one)
                   # Colors to be used in the gradient
                   col.min = 0,                       # Minimum scaled average expression threshold (everything smaller will be set to this)
                   col.max = 2.5,                        # Maximum scaled average expression threshold (everything larger will be set to this)
                   dot.min = 0,                          # The fraction of cells at which to draw the smallest dot (default is 0)
                   dot.scale = 6,                        # Scale the size of the points
                   group.by = NULL,              # How the cells are going to be grouped
                   split.by = NULL,                      # Whether to split the data (if you fo this make sure you have a different color for each variable)
                   scale = TRUE,                         # Whether the data is scaled
                   scale.by = "radius",                  # Scale the size of the points by 'size' or 'radius'
                   scale.min = NA,                       # Set lower limit for scaling
                   scale.max = NA )+                       # Set upper limit for scaling
  labs(x = NULL,                              # x-axis label
       y = NULL)+
  scale_color_viridis_c(option="F",begin=.4,end=0.9, direction = -1)+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.6)+
  #theme_linedraw()+
  guides(x =  guide_axis(angle = 90))+
  theme(axis.text.x = element_text(size = 8 , face = "italic"))+
  theme(axis.text.y = element_text(size = 9))+
  theme(legend.title = element_text(size = 9))+
  theme(legend.text = element_text(size = 8))+ 
  scale_y_discrete(limits = c("Ciliated 2",
                              "Ciliated 1",
                              "Selenop+/Gstm2high Secretory",
                              "Spdef+ Secretory",
                              "Fibroblast-like",
                              "Pax8low/Prom1+ Cilia-forming", 
                              "Slc1a3med/Sox9+ Cilia-forming",
                              "Cebpdhigh/Foxj1- Progenitor",
                              "Slc1a3+ Stem/Progenitor"))

ggsave(filename = "FIGs4_stem_dp.pdf", plot = stem_dp, width = 12, height = 6, dpi = 600)


x <- stem_dp$data

write.csv( x , 'stem_dp_data.csv')

#### Figure Supp 6: HGSC Driver Gene by Pseudotime ####

## Calculate Pseudotime Values ##

pseudo <- pseudotime(cds)

Distal_PHATE@meta.data$Pseudotime <- pseudo # Add to Seurat Metadata

## Subset Seurat Object ##

color_cells <- DimPlot(Distal_PHATE , reduction = "phate", 
                       cols = c("#B20224", #1
                                "#35EFEF", #2
                                "#00A1C6", #3
                                "#A374B5", #4
                                "#9000C6", #5
                                "#EA68E1", #6
                                "lightgrey", #7
                                "#2188F7", #8
                                "#F28D86"),
                       pt.size = 0.7,
                       shuffle = TRUE,
                       seed = 0,
                       label = FALSE)


## Psuedotime and Lineage Assignment ##

cellID <- rownames(Distal_PHATE@reductions$phate@cell.embeddings)
phate_embeddings <- Distal_PHATE@reductions$phate@cell.embeddings
pseudotime_vals <- Distal_PHATE@meta.data$Pseudotime

combined_data <- data.frame(cellID, phate_embeddings, pseudotime_vals)

# Calculate the Average PHATE_1 Value for Pseudotime Points = 0 #
avg_phate_1 <- mean(phate_embeddings[pseudotime_vals == 0, 1])

# Pseudotime Values lower than avge PHATE_1 Embedding will be Negative to split lineages
combined_data$Split_Pseudo <- ifelse(phate_embeddings[, 1] < avg_phate_1, -pseudotime_vals, pseudotime_vals)

# Define Lineage #
combined_data$lineage <- ifelse(combined_data$PHATE_1 < avg_phate_1, "Secretory",
                                ifelse(combined_data$PHATE_1 > avg_phate_1, "Ciliogenic", "Progenitor"))


Distal_PHATE$Pseudotime_Adj <- combined_data$Split_Pseudo
Distal_PHATE$Lineage <- combined_data$lineage

# Subset #

Pseudotime_Lineage <- subset(Distal_PHATE, 
                             idents = c("Secretory 1",
                                        "Secretory 2",
                                        "Msln+ Progenitor",
                                        "Slc1a3+/Sox9+ Cilia-forming",
                                        "Pax8+/Prom1+ Cilia-forming",
                                        "Progenitor",
                                        "Ciliated 1",
                                        "Ciliated 2"))


## Set Bins ##

bins <- cut_number(Pseudotime_Lineage@meta.data$Pseudotime_Adj , 40) # Evenly distribute bins 

Pseudotime_Lineage@meta.data$Bin <- bins # Metadata for Bins

## Set Idents to PSeudoime Bin ##

time_ident <- SetIdent(Pseudotime_Lineage, value = Pseudotime_Lineage@meta.data$Bin)

av.exp <- AverageExpression(time_ident, return.seurat = T)$RNA # Calculate Avg log normalized expression

# Calculates Average Expression for Each Bin #
# if you set return.seurat=T, NormalizeData is called which by default performs log-normalization #
# Reported as avg log normalized expression #


## Pseudotime Scale Bar ##

list <- 1:40
colors = c(rev(rainbow20),rainbow20)
df <- data.frame(data = list, color = colors)

pseudo_bar <- ggplot(df, aes(x = 1:40, y = 1, fill = color)) + 
  geom_bar(stat = "identity",position = "fill", color = "black", size = 0, width = 1) +
  scale_fill_identity() +
  theme_void()+ 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

ggsave(filename = "pseudo_bar.pdf", plot = pseudo_bar, width = 0.98, height = 0.19, dpi = 600)


## Plot HGSC driver gene list across pseudotime bin ##

features <- c("Trp53", "Brca1", "Brca2",	"Csmd3",	"Nf1",	"Fat3",	"Gabra6", "Rb1", "Apc",	"Lrp1b",
              "Prim2",	"Cdkn2a", "Crebbp",	"Wwox", "Ankrd11",	
              "Map2k4",	"Fancm",	"Fancd2",	"Rad51c",  "Pten")

# Create Bin List and expression of features #

bin_list <- unique(Pseudotime_Lineage@meta.data$Bin) 

plot_info <- as.data.frame(av.exp[features,]) # Call Avg Expression for features


z_score <- transform(plot_info, SD=apply(plot_info,1, mean, na.rm = TRUE))
z_score <- transform(z_score, MEAN=apply(plot_info,1, sd, na.rm = TRUE))

z_score1 <- (plot_info-z_score$MEAN)/z_score$SD



plot_info$y <- rownames(plot_info) # y values as features
z_score1$y <- rownames(plot_info)


plot_info <- gather(data = plot_info, x, expression, bin_list) #set plot
z_score1 <- gather(data = z_score1, x, z_score, bin_list) #set plot


# Create Cell Clusters DF #

Labeled_Pseudotime_Lineage <- RenameIdents(Pseudotime_Lineage, 
                                           'Secretory 1' = "Spdef+ Secretory", 
                                           'Progenitor' = "Slc1a3+ Stem/Progenitor", 
                                           'Msln+ Progenitor' = "Cebpdhigh/Foxj1- Progenitor",
                                           'Ciliated 1' = "Ciliated 1", 
                                           'Ciliated 2' = "Ciliated 2", 
                                           'Pax8+/Prom1+ Cilia-forming' = "Pax8low/Prom1+ Cilia-forming", 
                                           'Fibroblast-like' = "Fibroblast-like", #removed
                                           'Slc1a3+/Sox9+ Cilia-forming' = "Slc1a3med/Sox9+ Cilia-forming",
                                           'Secretory 2' = "Selenop+/Gstm2high Secretory")

cluster_table <- table(Labeled_Pseudotime_Lineage@active.ident, 
                       Labeled_Pseudotime_Lineage@meta.data$Bin)

clusters <- data.frame(cluster_table)

clusters <- clusters %>% 
  group_by(Var2) %>%
  mutate(Perc = Freq / sum(Freq))


# Create Pseudotime DF #

pseudotime_table <- table(seq(1, length(bin_list), 1), 
                          unique(Labeled_Pseudotime_Lineage@meta.data$Bin),
                          seq(1, length(bin_list), 1))

pseudotime_bins <- data.frame(pseudotime_table)  


# calculate max and min z-scores
max_z <- max(z_score1$z_score, na.rm = TRUE)
min_z <- min(z_score1$z_score, na.rm = TRUE)

# set color for outliers
outlier_color <- ifelse(z_score1$z_score > max_z | z_score1$z_score < min_z, ifelse(z_score1$z_score > 0, "#AD1F24", "#51A6DC"), "#e2e2e2")


## Plot Gene Expression ##

# Set different na.value options for positive and negative values
na_color_pos <- "#AD1F24" # color for positive NA values
na_color_neg <- "#51A6DC" # color for negative NA values

custom_bin_names <- c(paste0("S", 20:1), paste0("C", 1:20))

figure <- ggplot(z_score1, aes(x, y, fill = z_score)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradientn(colors=c("#1984c5", "#e2e2e2", "#c23728"), 
                       name = "Average Expression \nZ-Score", limits = c(-3,3), 
                       na.value = ifelse(is.na(z_score1) & z_score1 > 0, na_color_pos, 
                                         ifelse(is.na(z_score1) & z_score1 < 0, na_color_neg, "grey50")),
                       oob = scales::squish)+
  scale_x_discrete(limits= sort(bin_list) , labels= custom_bin_names)+
  scale_y_discrete(limits= rev(features))+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),                                      # Text size throughout the plot
        axis.text.x = element_text(color = 'black', angle = 0, hjust = 0.5, size = 10, face = "bold"),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold.italic"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-0.5,1,1,1), "cm"))


## Plot Cluster Percentage ##


`Spdef+ Secretory` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Spdef+ Secretory")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(1,1,1,1), "cm"))

`Selenop+/Gstm2high Secretory` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Selenop+/Gstm2high Secretory")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

`Cebpdhigh/Foxj1- Progenitor` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Cebpdhigh/Foxj1- Progenitor")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

`Slc1a3+ Stem/Progenitor` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Slc1a3+ Stem/Progenitor")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

`Slc1a3med/Sox9+ Cilia-forming` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Slc1a3med/Sox9+ Cilia-forming")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

`Pax8low/Prom1+ Cilia-forming` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Pax8low/Prom1+ Cilia-forming")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

`Ciliated 1` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Ciliated 1")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))

`Ciliated 2` <- ggplot(clusters, aes(Var2, Var1, fill = Perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1) +
  scale_fill_gradient2(low="white", high="#000000", mid = "white", midpoint = 0, name = "Percentage")+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Ciliated 2")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust = 1, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))


## Plot Pseudotime Color ##

list <- 1:40
colors = c(rev(rainbow20),rainbow20)
df <- data.frame(data = list, color = colors)


binning <- ggplot(df, aes(x = 1:40, y = 1, fill = color)) + 
  geom_bar(stat = "identity",position = "fill", color = "black", size = 1, width = 1) +
  scale_fill_identity() +
  theme_void()+ 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())+
  scale_x_discrete(limits= sort(bin_list) , labels= seq(1, length(bin_list), 1))+
  scale_y_discrete(limits= "Pseudotime Bin ")+
  theme(panel.background = element_blank())+
  labs(title = "Expression of Genes by Pseudotime Bin" ,
       x = element_blank(),
       y = element_blank())+
  theme(text = element_text(size = 12, face = "bold"),# Text size throughout the plot
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),    # Text color, angle, and horizontal adjustment on x-axis 
        axis.text.y = element_text(color = 'black', hjust =1, vjust = .75, size = 14, face = "bold"))+
  theme(plot.title = element_blank(),
        plot.margin=unit(c(-1.25,1,1,1), "cm"))


### Combine Plots ###


psuedotime_lineage <- ggarrange(`Spdef+ Secretory`,
                                `Selenop+/Gstm2high Secretory`,
                                `Cebpdhigh/Foxj1- Progenitor`,
                                `Slc1a3+ Stem/Progenitor`,
                                `Slc1a3med/Sox9+ Cilia-forming`,
                                `Pax8low/Prom1+ Cilia-forming`,
                                `Ciliated 1`,
                                `Ciliated 2`,
                                `binning`,
                                figure , ncol=1,
                                heights = c(2, 2, 2, 2, 2, 2, 2, 2, 2, (2*length(features)),
                                            widths = c(3)),
                                padding = unit(0.01))

ggsave(filename = "FIGs6_psuedotime_driver_gene.pdf", plot = psuedotime_lineage, width = 18, height = 9, dpi = 600)


write.csv(z_score1 , 'cancer_pseudotime.csv')

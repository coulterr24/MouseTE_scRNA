#### Figure Supp 3: Doublet detection of fibroblast and epithelial markers ####

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


#### Load Distal Epithelial Dataset ####

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

table(Epi_Named@active.ident)


#### Plot Compilation ####

feature_scatter <- FeatureScatter( Fibroblast, "Krt8","Col1a1",
                                   cols = c("#35EFEF", #1
                                            "#00A1C6", #2
                                            "#2188F7", #3
                                            "#EA68E1", #4
                                            "#59D1AF", #5
                                            "#B20224", #6
                                            "#F28D86", #7
                                            "#A374B5", #8
                                            "#9000C6"))

x <- DotPlot(Epi_Named , features = c("Krt8" , "Col1a1"))
write.csv(x$data , "doublet_data.csv")



Fibroblast <- subset(Epi_Named, 
                     idents = c("Fibroblast-like"))



custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( Fibroblast, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))


ggsave(filename = "20240221_Fibroblast_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)

## Stem ##

c("Cebpdhigh/Foxj1- Progenitor",
  "Slc1a3med/Sox9+ Cilia-forming",
  "Pax8low/Prom1+ Cilia-forming",
  "Spdef+ Secretory",
  "Selenop+/Gstm2high Secretory",
  "Ciliated 1",
  "Ciliated 2")


Stem <- subset(Epi_Named, 
               idents = c("Slc1a3+ Stem/Progenitor"))


custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( Stem, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))


ggsave(filename = "20240221_Slc1a3_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)



## Prog ##

c("Slc1a3med/Sox9+ Cilia-forming",
  "Pax8low/Prom1+ Cilia-forming",
  "Spdef+ Secretory",
  "Selenop+/Gstm2high Secretory",
  "Ciliated 1",
  "Ciliated 2")


Prog <- subset(Epi_Named, 
               idents = c("Cebpdhigh/Foxj1- Progenitor"))


custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( Prog, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_Cebpd_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)



## Cilia-forming ##

c("Pax8low/Prom1+ Cilia-forming",
  "Spdef+ Secretory",
  "Selenop+/Gstm2high Secretory",
  "Ciliated 1",
  "Ciliated 2")


trans <- subset(Epi_Named, 
                idents = c("Slc1a3med/Sox9+ Cilia-forming"))

custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( trans, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_SlcCilia_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)


## Pax8 Cilia-forming ##

c("Spdef+ Secretory",
  "Selenop+/Gstm2high Secretory",
  "Ciliated 1",
  "Ciliated 2")


cancer <- subset(Epi_Named, 
                 idents = c("Pax8low/Prom1+ Cilia-forming"))

custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( cancer, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_Prom1_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)



## Spdef Secretory ##

c("Selenop+/Gstm2high Secretory",
  "Ciliated 1",
  "Ciliated 2")


sec1 <- subset(Epi_Named, 
               idents = c("Spdef+ Secretory"))

custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( sec1, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_Spdef_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)


## Selenop Secretory ##

c("Ciliated 1",
  "Ciliated 2")


sec2 <- subset(Epi_Named, 
               idents = c("Selenop+/Gstm2high Secretory"))

custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( sec2, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_Selenop_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)


## Ciliated 1 ##

c("Ciliated 2")


cil1 <- subset(Epi_Named, 
               idents = c("Ciliated 1"))

custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( cil1, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_Ciliated_1_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)



## Ciliated 2 ##


cil2 <- subset(Epi_Named, 
               idents = c("Ciliated 2"))

custom_labels <- function(x) {
  ifelse(x %% 1 == 0, as.character(x), "")
}

feature_scatter <- FeatureScatter( cil2, "Krt8","Col1a1",
                                   cols = "black")+  # Scatter plot
  NoLegend()+
  labs(title = NULL)+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5),
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(color = 'black', size = 12),
        axis.title.x = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'black'),
        axis.text.y = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'black'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels) +  # Set breaks every 0.5 units on x-axis
  scale_y_continuous(limits = c(0,5) , breaks = seq(0, 10, by = 0.5),
                     labels = custom_labels)   # Set breaks every 0.5 units on y-axis



plot_data <- as.data.frame(feature_scatter$data)


# Create density plot for x-axis
density_x <- ggplot(plot_data, aes(x = Krt8 , fill = 'black')) +
  geom_density() +
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color = 'white', size = 12),
        axis.title.y = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.y = element_line(color = 'white'),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  scale_y_continuous(labels = function(y) sprintf("%.0f", y))+ 
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,4) , breaks = seq(0, 10, by = 0.5))  # Set breaks every 0.5 units on x-axis



# Create density plot for y-axis
density_y <- ggplot(plot_data, aes(x = Col1a1 , fill = 'black')) +
  geom_density() +
  coord_flip() +  # Flip axes for y-density plot
  theme(axis.line = element_line(color='white'),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(color = 'white', size = 12),
        axis.title.x = element_text(color = 'white', size = 14, face = "bold.italic"),
        axis.ticks.x = element_line(color = 'white'),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "mm"))+
  NoLegend()+
  scale_fill_manual(values = c("grey"))+
  scale_x_continuous(limits = c(0,5))+
  scale_y_continuous(limits = c(0,1) , breaks = seq(0, 10, by = 0.5))   # Set breaks every 0.5 units on y-axis

# Arrange plots

top_row <- cowplot::plot_grid(
  density_x, NULL,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(0.5,0.5))

bottom_row <- cowplot::plot_grid(
  feature_scatter,density_y,
  nrow = 1, rel_widths = c(3, 1), rel_heights = c(3,3))

combined_plot <- cowplot::plot_grid(
  top_row , bottom_row,
  nrow = 2 , rel_widths = c(3, 1), rel_heights = c(0.25,3))

ggsave(filename = "20240221_Ciliated_2_Doublet.pdf", plot = combined_plot, width = 12, height = 12, dpi = 600)






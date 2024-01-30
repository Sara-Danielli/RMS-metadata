library(dplyr)
library(Seurat)
library(ggplot2)

## color palette
col_cluster_names_aggregate <- c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF' , '#D8D8E0FF', '#B497E7FF')
names(col_cluster_names_aggregate) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis', 'IFN')

col_cluster_names_aggregate_integrated <- c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF", "#8E0152FF", '#D8D8E0FF')
names(col_cluster_names_aggregate_integrated) <- c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis')

## Style plots

# Barplot function
theme_ggplot = theme(legend.position = "none",
                     plot.title = element_text(hjust=0.5, face="bold"),
                     panel.border = element_blank(),
                     plot.background = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
                     axis.text.y = element_text(size=12, colour="black"),
                     axis.title=element_text(size=12),
                     strip.text = element_text(size = 13, face = "bold"))

theme_ggplot_legend = theme(plot.title = element_text(hjust=0.5, face="bold"),
                            panel.border = element_blank(),
                            plot.background = element_blank(),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size=12, angle = 45, vjust = 1, hjust=1, colour="black"),
                            axis.text.y = element_text(size=12, colour="black"),
                            axis.title=element_text(size=12),
                            strip.text = element_text(size = 13, face = "bold"))

## Style plots

theme_vln <- theme(panel.border = element_blank(),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
                   axis.text.y = element_text(size=14, colour="black"),
                   axis.title=element_text(size=14),
                   plot.title = element_text(size=14, face="bold")) 


theme_cellstate_plot <-  theme(panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               plot.background = element_blank(),
                               panel.background = element_blank(),
                               axis.line = element_line(colour = "black"),
                               axis.title=element_text(size=16),
                               axis.text.x = element_text(size=16, vjust = 0.5, colour="black"),
                               axis.text.y = element_text(size=16, colour="black"),
                               plot.title = element_text(hjust = 0.5, angle = 0, size = 16, face = "bold", vjust = 1),
                               legend.text = element_text(size = 12),
                               legend.title = element_text(size = 15, face = "bold"),
                               legend.box.background = element_rect(color="black", size=1),
                               strip.text = element_text(size = 15, face = "bold"))

# Barplot function
plot_bar <- function(seurat_obj, x_var, y_var, colors){
  ggplot(seurat_obj@meta.data, aes(x_var, fill = y_var)) +
    scale_fill_manual(values = colors) + 
    geom_bar(position = "fill", color="black") +
    labs (y='Proportion', x='') + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=14, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=14, colour="black"),
          axis.title=element_text(size=14),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) 
}

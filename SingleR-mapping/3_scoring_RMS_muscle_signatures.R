rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(SingleCellExperiment)
library(readxl)
library(data.table)
library(paletteer)
library(RColorBrewer)
#library(scater)
#library(scRNAseq)
library(SingleR)

# Set up environment ------------------------------------------------
setwd('/Volumes/Sara_PhD/scRNAseq_data/')
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
base_dir <- file.path('/n/scratch/users/s/sad167/RMS')

# Load data ------------------------------------------------
# load data   
PDX.integrated <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))
muscle_signatures <- readRDS(file.path(base_dir, 'data/Xi/7_markers.rds'))
  
# Score tumor programs for gene expression programs identified in Xi et al. --------------------
# Score RMS for muscle signatures
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = muscle_signatures, name = names(muscle_signatures))
PDX.integrated <- ScaleData(PDX.integrated)
# rename metadata names of scores
# identify number of first column with metadata scores
col_start <- length(colnames(PDX.integrated@meta.data)) - length(names(muscle_signatures)) +1
# identify number of last column with metadata scores
col_end <- length(colnames(PDX.integrated@meta.data))
# rename columns with score name
colnames(PDX.integrated@meta.data)[col_start:col_end] <- names(muscle_signatures)

# Save file
saveRDS(PDX.integrated, file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_20230202_scoring100_20240115.rds"))

 
# Plot results  --------------------
DotPlot(PDX.integrated, features = names(muscle_signatures), 
        assay = 'RNA', 
        col.min = 0, 
        cols = c("white", "red3"),
        scale = FALSE)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/5_DotPlot_scores.pdf", 
       width=5, height=4, dpi=300)
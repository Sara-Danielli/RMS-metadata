rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(data.table)
library(paletteer)
library(SeuratObject)
library(tidyverse)
library(cowplot)
library(ggpubr)

base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
data_dir <- file.path(base_dir, 'RMS/Patel_2022/FFPE_RNAseq')
plot_dir <- file.path(base_dir, 'output/metadata/Pseudobulk')

if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}


##################################################################
#  (0) Loading Patel FFPE RNAseq datasets:  #
##################################################################

## ! Excel matrix contained some numbers that were not integers --> approximated

# load TPM data
dataset <- read.delim(file.path(data_dir, "RMS13_FFPE_RNA-seq_tpm_data.txt"))

rownames(dataset) <- make.names(dataset$geneSymbol, unique = TRUE)

# set gene names as row header
dataset2 <- dataset[, -1]
dataset3 <- dataset2[, -1]
dataset4 <- dataset3[, -1]
dataset5 <- dataset4[, -1]

# load metadata 
FFPE_meta <- read_excel(file.path(data_dir,"metadata.xlsx"), col_names=TRUE)
  # set first column as row header
FFPE_meta2 <- FFPE_meta %>% remove_rownames %>% column_to_rownames(var="patient_id")

# Create Seurat object
FFPE_seurat <- CreateSeuratObject(counts = dataset5)
FFPE_seurat <- AddMetaData(FFPE_seurat, metadata = FFPE_meta2)

# log-transform (as data are already TPM, do not perform "NormalizeData")
counts <- GetAssayData(FFPE_seurat, slot = "counts")
counts_log <-log1p(counts)

# Re-create Seurat object with log-transformed counts
FFPE_seurat <- CreateSeuratObject(counts = counts_log)
FFPE_seurat <- AddMetaData(FFPE_seurat, metadata = FFPE_meta2)

# Reoder levels
FFPE_seurat$treatment <- factor(x = FFPE_seurat$treatment, levels = c('diagnostic biopsy', 'delayed resection'))

##################################################################
#  (1) Score datasets with gene signatures:  #
##################################################################

# Load common signatures
  markers <- list()
  markers$Common_differentiated <- read_excel("/Volumes/Sara_PhD/scRNAseq_data/lists/Common_differentiated.xlsx", col_names=FALSE)
  markers$Common_proliferative <- read_excel("/Volumes/Sara_PhD/scRNAseq_data/lists/Common_proliferative.xlsx", col_names=FALSE)
  markers$Common_Stemcell <- read_excel("/Volumes/Sara_PhD/scRNAseq_data/lists/Common_Stemcell.xlsx", col_names=FALSE)
  

# Add scores using "Addmodulescore" fct
  FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = markers$Common_differentiated, ctrl = 100, name = "Common_differentiated")
  FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = markers$Common_proliferative, ctrl = 100, name = "Common_proliferative")
  FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = markers$Common_Stemcell, ctrl = 100, name = "Common_Stemcell")
 
  FFPE_seurat$DS_score_Addmodule <- FFPE_seurat$Common_differentiated1 - FFPE_seurat$Common_Stemcell1
  

# Scale data
ScaleData(FFPE_seurat)

      
### Vln plot Addmodulescore
      compare_means(Common_proliferative1 ~ treatment,  data = FFPE_seurat[[]], paired = TRUE)
      p2 <- VlnPlot(FFPE_seurat,
              features = 'Common_proliferative1',
              group.by = 'treatment',
              pt.size=0,
              cols = c('lightseagreen', 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype,group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(paired = TRUE, size=2) + NoLegend()
      ggsave(file.path(plot_dir, "3_Vln_plot_prolifer_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      
      
      compare_means(Common_differentiated1 ~ treatment,  data = FFPE_seurat[[]], paired = TRUE)
      p3 <- VlnPlot(FFPE_seurat,
              features = 'Common_differentiated1',
              group.by = 'treatment',
              pt.size=0,
              cols = c('lightseagreen', 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype, group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(paired = TRUE, size=2) + NoLegend()
      ggsave(file.path(plot_dir, "2_Vln_plot_differentiated_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      
      
      
      compare_means(Common_Stemcell1 ~ treatment,  data = FFPE_seurat[[]])
      p1 <- VlnPlot(FFPE_seurat,
              features = 'Common_Stemcell1',
              group.by = 'treatment',
              pt.size=0,
              cols = c('lightseagreen', 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype, group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(paired = TRUE, size=2) + NoLegend()
      ggsave(file.path(plot_dir, "2_Vln_plot_progenitor_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      
      
      compare_means(DS_score_Addmodule ~ treatment,  data = FFPE_seurat[[]])
      VlnPlot(FFPE_seurat,
              features = 'DS_score_Addmodule',
              group.by = 'treatment',
              pt.size=0,
              cols = c('lightseagreen', 'gray96', 'gray96', 'maroon'),
      ) + 
        geom_point(aes(fill=FFPE_seurat@meta.data$subtype, group=FFPE_seurat@meta.data$orig.ident),size=3,shape=21, position = position_dodge(0.2)) +
        geom_line(aes(group = FFPE_seurat@meta.data$orig.ident), position = position_dodge(0.2)) +
        stat_compare_means(paired = TRUE, size=2) + NoLegend()
      ggsave(file.path(plot_dir, "4_Vln_plot_muscle-score_treatment.pdf"), width=2.5, height=4, dpi=300) 
      
      plot_grid(p1, p2, p3, ncol = 3, align = "h", axis = "tb")
      ggsave(file.path(plot_dir,"0_VlnPlot.pdf"), width=6, height=3.5, dpi=300)
      
 
      
           5## Export data tables
md <- FFPE_seurat@meta.data %>% as.data.table
table_PDX <- md[, .N, by = c("name", "subtype", "Common_differentiated1", "Common_proliferative1", "Common_Stemcell1", "DS_score_Addmodule")] 
write.csv(table_PDX, file.path(plot_dir, "Signatures.csv"))
      
 
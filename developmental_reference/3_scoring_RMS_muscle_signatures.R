rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratObject)
library(readxl)
library(data.table)
library(paletteer)
library(RColorBrewer)

# Set up environment ------------------------------------------------
#base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
base_dir <- file.path('/n/scratch/users/s/sad167/RMS')

analysis_dir <- file.path(base_dir, 'analysis')
plot_dir <- file.path(base_dir, 'analysis/plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
write_dir <- file.path(base_dir, 'analysis/data')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}

# color palette
col_cluster_names_aggregate <- c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF' , '#D8D8E0FF', '#B497E7FF')
names(col_cluster_names_aggregate) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis', 'IFN')

# Load data ------------------------------------------------
#PDX.integrated <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))
P3F1 <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_20230713.rds"))
P7F1 <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_20230713.rds"))
ERMS <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_ERMS_20230713.rds"))

# rename P3F1 (problem with current IDs)
Idents(P3F1) = "cluster_names"
new.cluster.ids.aggregate <- c('Progenitor', 'Proliferative', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')
names(new.cluster.ids.aggregate) <- levels(P3F1)
P3F1<- RenameIdents(P3F1, new.cluster.ids.aggregate)
P3F1[["cluster_names_aggregate"]] <- Idents(object = P3F1)
levels(P3F1) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')



seurat_objects <- list(#PDX.integrated, 
  P3F1, P7F1, ERMS)
names(seurat_objects) <- c('P3F1', 'P7F1', 'FNRMS')

muscle_signatures <- readRDS(file.path(base_dir, 'data/Xi/7_markers.rds'))

# Score tumor programs for gene expression programs identified in Xi et al. --------------------
seurat_objects <- lapply(seurat_objects, FUN = function(x) {
  # Score RMS for muscle signatures
  x <- AddModuleScore(object = x, assay = 'RNA', features = muscle_signatures, name = names(muscle_signatures))
  x <- ScaleData(x)
})

# rename metadata names of scores
col_start <- list()
col_end <- list()
for (i in seq_along(seurat_objects)){
  # identify number of first column with metadata scores
  col_start[[i]] <- length(colnames(seurat_objects[[i]]@meta.data)) - length(names(muscle_signatures)) +1
  # identify number of last column with metadata scores
  col_end[[i]] <- length(colnames(seurat_objects[[i]]@meta.data))
  # rename columns with score name
  colnames(seurat_objects[[i]]@meta.data)[col_start[[i]]:col_end[[i]]] <- names(muscle_signatures)
}

# Plot results  --------------------
for (i in seq_along(seurat_objects)){
  Idents(seurat_objects[[i]]) <- "cluster_names_aggregate"
  
  # UMAP plot
  DimPlot(seurat_objects[[i]], reduction = "umap_rpca", cols = col_cluster_names_aggregate, pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 
  ggsave(file.path(plot_dir, paste0("1_UMAP_", names(seurat_objects)[i], ".pdf")), width=5.5, height=4.5, dpi=300)
  
  # Dotplot scores
  DotPlot(seurat_objects[[i]], 
          features = names(muscle_signatures), 
          #split.by = 'subtype',
          assay = 'RNA', 
          cols = c("white", "red3"),
          col.min = 0, 
          scale = FALSE)  + 
    #scale_colour_distiller(palette="RdBu") +
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=12, colour="black"),
          axis.title=element_blank(),
          legend.text = element_text(size = 12),
          #legend.title = element_blank()
    ) 
  ggsave(file.path(plot_dir, paste0("2_DotPlot_scores_", names(seurat_objects)[i], ".pdf")), width=5.5, height=4.5, dpi=300)
}


#  ------------------------------------------------ ------------------------------------------------
# Repeat same with combined FPRMS and FNRMS atlas  ------------------------------------------------
#  ------------------------------------------------ ------------------------------------------------

# Color palette
col_cluster_names_aggregate <- c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF", "#8E0152FF", '#D8D8E0FF')
names(col_cluster_names_aggregate) <- c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis')

# Load data ------------------------------------------------
PDX.integrated <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))

muscle_signatures <- readRDS(file.path(base_dir, 'data/Xi/7_markers.rds'))
seurat_objects <- list(PDX.integrated)
names(seurat_objects) <- c('PDX_integrated')

# Score tumor programs for gene expression programs identified in Xi et al. --------------------
seurat_objects <- lapply(seurat_objects, FUN = function(x) {
  # Score RMS for muscle signatures
  x <- AddModuleScore(object = x, assay = 'RNA', features = muscle_signatures, name = names(muscle_signatures))
  x <- ScaleData(x)
})

# rename metadata names of scores
col_start <- list()
col_end <- list()
for (i in seq_along(seurat_objects)){
  # identify number of first column with metadata scores
  col_start[[i]] <- length(colnames(seurat_objects[[i]]@meta.data)) - length(names(muscle_signatures)) +1
  # identify number of last column with metadata scores
  col_end[[i]] <- length(colnames(seurat_objects[[i]]@meta.data))
  # rename columns with score name
  colnames(seurat_objects[[i]]@meta.data)[col_start[[i]]:col_end[[i]]] <- names(muscle_signatures)
  
  # Save file
  #saveRDS(seurat_objects[obj], file.path(write_dir, "Danielli_Patel_Langenau_RPCA_20230202_scoring100_20240115.rds"))
}

# Plot results  --------------------
for (i in seq_along(seurat_objects)){
  Idents(seurat_objects[[i]]) <- "cluster_names_aggregate"
  
  # UMAP plot
  DimPlot(seurat_objects[[i]], reduction = "umap_rpca", cols = col_cluster_names_aggregate, pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 
  ggsave(file.path(plot_dir, paste0("1_UMAP_", names(seurat_objects)[i], ".pdf")), width=5.5, height=4.5, dpi=300)
  
  # Dotplot scores
  DotPlot(seurat_objects[[i]], 
          features = names(muscle_signatures), 
          #split.by = 'subtype',
          assay = 'RNA', 
          cols = c("white", "red3"),
          scale = FALSE,
          col.min = 0)  + 
    #scale_colour_distiller(palette="RdBu") +
    theme(axis.line = element_line(colour = "black"),
          axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
          axis.text.y = element_text(size=12, colour="black"),
          axis.title=element_blank(),
          legend.text = element_text(size = 12),
          #legend.title = element_blank()
    ) 
  ggsave(file.path(plot_dir, paste0("2_DotPlot_scores_", names(seurat_objects)[i], ".pdf")), width=5.5, height=4.5, dpi=300)
}



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
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
resource_dir <- file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources')
source(file.path(resource_dir, 'Plot_style_v2.R'))

# Load data FP-RMS (PAX3::FOXO1) ------------------------------------------------
P3F1 <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_20230713.rds"))

# rename P3F1 (problem with current IDs)
Idents(P3F1) = "cluster_names"
new.cluster.ids.aggregate <- c('Progenitor', 'Proliferative', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')
names(new.cluster.ids.aggregate) <- levels(P3F1)
P3F1<- RenameIdents(P3F1, new.cluster.ids.aggregate)
P3F1[["cluster_names_aggregate"]] <- Idents(object = P3F1)
levels(P3F1) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')

# UMAP plot to confirm correct assignment
DimPlot(P3F1, reduction = "umap_rpca", cols = col_cluster_names_aggregate,  pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 

# clean up metadata
P3F1[['orig.ident']] <- NULL
P3F1[['old.ident']] <- NULL
P3F1[['Wei_FN_genes1']] <- NULL
P3F1[['Wei_FP_genes1']] <- NULL
P3F1[['integrated_snn_res.0.3']] <- NULL
P3F1[['cluster_names']] <- NULL

colnames(P3F1@meta.data)

colnames(P3F1@meta.data)[23] <- 'Differentiated score'
colnames(P3F1@meta.data)[24] <- 'Progenitor score'
colnames(P3F1@meta.data)[25] <- 'Proliferative score'
colnames(P3F1@meta.data)[26] <- 'Muscle lineage score'
colnames(P3F1@meta.data)[28] <- 'Cluster assignment'

head(P3F1@meta.data)

# save dataset
saveRDS(P3F1, file.path(base_dir, 'write/FPRMS_PAX3FOXO1_final_20240130.rds'))

rm(P3F1)


# Load data FP-RMS (PAX7::FOXO1) ------------------------------------------------
P7F1 <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_20230713.rds"))

# UMAP plot to confirm correct assignment
DimPlot(P7F1, reduction = "umap_rpca", group.by = 'cluster_names_aggregate', cols = col_cluster_names_aggregate,  pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 

# clean up metadata
P7F1[['orig.ident']] <- NULL
P7F1[['old.ident']] <- NULL
P7F1[['Wei_FN_genes1']] <- NULL
P7F1[['Wei_FP_genes1']] <- NULL
P7F1[['integrated_snn_res.0.2']] <- NULL
P7F1[['cluster_names']] <- NULL

colnames(P7F1@meta.data)

colnames(P7F1@meta.data)[23] <- 'Differentiated score'
colnames(P7F1@meta.data)[24] <- 'Progenitor score'
colnames(P7F1@meta.data)[25] <- 'Proliferative score'
colnames(P7F1@meta.data)[26] <- 'Muscle lineage score'
colnames(P7F1@meta.data)[28] <- 'Cluster assignment'

head(P7F1@meta.data)

# save dataset
saveRDS(P7F1, file.path(base_dir, 'write/FPRMS_PAX7FOXO1_final_20240130.rds'))

rm(P7F1)

# Load data FN-RMS ------------------------------------------------
FNRMS <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ERMS_20230713.rds"))

# UMAP plot to confirm correct assignment
DimPlot(FNRMS, reduction = "umap_rpca", group.by = 'cluster_names_aggregate', cols = col_cluster_names_aggregate,  pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 

# clean up metadata
FNRMS[['orig.ident']] <- NULL
FNRMS[['old.ident']] <- NULL
FNRMS[['Wei_FN_genes1']] <- NULL
FNRMS[['Wei_FP_genes1']] <- NULL
FNRMS[['integrated_snn_res.0.2']] <- NULL
FNRMS[['cluster_names']] <- NULL

colnames(FNRMS@meta.data)

colnames(FNRMS@meta.data)[23] <- 'Progenitor score'
colnames(FNRMS@meta.data)[24] <- 'Proliferative score'
colnames(FNRMS@meta.data)[25] <- 'Differentiated score'
colnames(FNRMS@meta.data)[26] <- 'Muscle lineage score'
colnames(FNRMS@meta.data)[28] <- 'Cluster assignment'

head(FNRMS@meta.data)

# save dataset
saveRDS(FNRMS, file.path(base_dir, 'write/FNRMS_final_20240130.rds'))

rm(FNRMS)



# Load data ------------------------------------------------
PDX.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))

# UMAP plot to confirm correct assignment
DimPlot(PDX.integrated, reduction = "umap_rpca", group.by = 'cluster_names_aggregate', cols = col_cluster_names_aggregate_integrated,  pt.size = 4, raster=TRUE, shuffle=TRUE, raster.dpi = c(1012, 1012)) 

# clean up metadata
PDX.integrated[['orig.ident']] <- NULL
PDX.integrated[['old.ident']] <- NULL
PDX.integrated[['Wei_mesenchymal1']] <- NULL
PDX.integrated[['Wei_proliferative1']] <- NULL
PDX.integrated[['Wei_muscle1']] <- NULL
PDX.integrated[['Wei_FN1']] <- NULL
PDX.integrated[['Wei_FP1']] <- NULL
PDX.integrated[['Danielli_MuSC1']] <- NULL
PDX.integrated[['Danielli_cycling1']] <- NULL
PDX.integrated[['Danielli_differentiated1']] <- NULL
PDX.integrated[['Patel_mesoderm1']] <- NULL
PDX.integrated[['Patel_myoblast1']] <- NULL
PDX.integrated[['Patel_myocyte1']] <- NULL
PDX.integrated[['integrated_snn_res.0.3']] <- NULL
PDX.integrated[['cluster_names']] <- NULL

colnames(PDX.integrated@meta.data)

colnames(PDX.integrated@meta.data)[23] <- 'Differentiated score'
colnames(PDX.integrated@meta.data)[24] <- 'Progenitor score'
colnames(PDX.integrated@meta.data)[25] <- 'Proliferative score'
colnames(PDX.integrated@meta.data)[26] <- 'Muscle lineage score'
colnames(PDX.integrated@meta.data)[28] <- 'Cluster assignment'

head(PDX.integrated@meta.data)

# save dataset
saveRDS(PDX.integrated, file.path(base_dir, 'write/RMS_atlas_final_20240130.rds'))

rm(PDX.integrated)

rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(writexl)
library(data.table)
library(viridis)
library(future)
library(paletteer)
library(SeuratObject)
library(future)
library(clustree)
plan("multicore", workers = 4)

base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
source(file.path(base_dir, "codes/MANUSCRIPT_INTEGRATION/metadata/FINAL/Functions.R"))

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/ERMSanalysis')
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## read file
ERMS.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ERMS_scores_20230215.rds"))
DefaultAssay(ERMS.integrated) <- "integrated"


## Order factors
ERMS.integrated$origin <- factor(x = ERMS.integrated$origin, levels = c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al."))
ERMS.integrated$model <- factor(x = ERMS.integrated$model, levels = c("Patient", "O-PDX", "Primary culture", "Cell line"))
ERMS.integrated$subtype <- factor(x = ERMS.integrated$subtype, levels = c("eRMS", "aRMS"))
ERMS.integrated$name <- factor(x = ERMS.integrated$name, levels = c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
                                                                  'aRMS-4',  'aRMS-5', 'KFR', 'Mast118', 
                                                                  'Mast95', 'MSK72117', 'MSK72117_SC', 
                                                                  'MSK82489', 'Rh4',  'Rh41',   'RMS', 
                                                                  'SJRHB010468_D1',  'SJRHB010468_X1', 'SJRHB013757_D2', 
                                                                  'SJRHB013757_X1', 'SJRHB013759_A1',  'SJRHB013759_A2', 
                                                                  'SJRHB013759_X14','SJRHB013759_X15', 'SJRHB031320_D1', 
                                                                  'SJRHB031320_X1', 'SJRHB046156_A1', 'SJRHB046156_X1', 
                                                                  '20696','21202', '29806', 
                                                                  'eRMS-1.1','eRMS-1.2',  'eRMS-2.1', 'eRMS-2.2', 
                                                                  'eRMS-3.2','eRMS-4',  'eRMS-8.1', 'eRMS-8.2', 
                                                                  'eRMS-8.3', 'Mast111','Mast139',  'Mast139_SC', 
                                                                  'Mast39',   'Mast85_r1','Mast85_r2', 
                                                                  'Mast85_r2_SC', 'MSK74711', 'RD', 'SJRHB000026_R2',  'SJRHB000026_R3',  'SJRHB000026_X1', 
                                                                  'SJRHB000026_X2', 'SJRHB010927_D1', 'SJRHB010927_X1', 
                                                                  'SJRHB010928_R1', 'SJRHB010928_X1',  'SJRHB011_D', 
                                                                  'SJRHB011_X', 'SJRHB012_R', 'SJRHB012_S', 'SJRHB012_Y', 
                                                                  'SJRHB012_Z', 'SJRHB012405_D1', 'SJRHB012405_X1',  'SJRHB013758_D1', 
                                                                  'SJRHB013758_D2', 'SJRHB013758_X1', 'SJRHB013758_X2',  'SJRHB030680_R1', 
                                                                  'SJRHB030680_X1', 'SJRHB049189_D1', 'SJRHB049189_X1'))




## Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
ERMS.integrated <- CellCycleScoring(ERMS.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#ERMS.integrated <- ScaleData(ERMS.integrated)

## Define high-cycling vs low-cycling cells
Cycling <- WhichCells(ERMS.integrated, expression = S.Score > 0 | G2M.Score > 0)
nonCycling <- WhichCells(ERMS.integrated, expression = S.Score <= 0 & G2M.Score <= 0)

# Add cycling info to metadata column called 'Cycling_prop'
ERMS.integrated$Cycling_prop <- ifelse(colnames(ERMS.integrated) %in% Cycling, "Cycling", "Non-cycling")


## Run PCA
ERMS.integrated <- RunPCA(ERMS.integrated, reduction.name = "pca_rpca", reduction.key = "pca_rpca")

## Elbow plot to define number of PCs to use for downstream analysis    
ElbowPlot(ERMS.integrated, reduction = 'pca_rpca')
ggsave(file.path(analysis_dir, "1_Elbowplot_rpca.pdf"), width=6, height=4, dpi=300)

## Run UMAP and find clusters
ERMS.integrated <- FindNeighbors(ERMS.integrated, reduction = "pca_rpca", dims = 1:13)

# Define optimal number of clusters with clustree
PDX_clusters <- FindClusters(ERMS.integrated, resolution = seq(0, 0.5, 0.1))
clustree(PDX_clusters)
ggsave(file.path(analysis_dir,"2_clustre.pdf"), width=6, height=8, dpi=300)

ERMS.integrated <- FindClusters(ERMS.integrated, resolution = 0.2)

ERMS.integrated <- RunUMAP(ERMS.integrated, dims = 1:13, reduction = "pca_rpca",  reduction.key = "umap_rpca", reduction.name = "umap_rpca")


## change names of clusters
levels(ERMS.integrated)
new.cluster.ids <- c('Ground', 'DNA replication', '2-Ground', 'Proliferative', '4-Progenitor', 'Differentiated', '6-Progenitor', 'IFN')
names(new.cluster.ids) <- levels(ERMS.integrated)
ERMS.integrated <- RenameIdents(ERMS.integrated, new.cluster.ids)
ERMS.integrated[["cluster_names"]] <- Idents(object = ERMS.integrated)
Idents(ERMS.integrated) = "cluster_names"
ERMS.integrated$cluster_names <- factor(x = ERMS.integrated$cluster_names, 
                                        levels = c('4-Progenitor', '6-Progenitor', 'Proliferative', 'Differentiated', 'IFN', 'DNA replication', 'Ground',  '2-Ground'))

Idents(ERMS.integrated) = 'cluster_names'
levels(ERMS.integrated)
new.cluster.ids.aggregate <- c('Progenitor', 'Progenitor', 'Proliferative', 'Differentiated', 'IFN', 'Proliferative', 'Ground', 'Ground')
names(new.cluster.ids.aggregate) <- levels(ERMS.integrated)
ERMS.integrated <- RenameIdents(ERMS.integrated, new.cluster.ids.aggregate)
ERMS.integrated[["cluster_names_aggregate"]] <- Idents(object = ERMS.integrated)
levels(ERMS.integrated) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'IFN')


Idents(ERMS.integrated) = ERMS.integrated$cluster_names

## Find markers
ERMS.integrated.markers <- FindAllMarkers(ERMS.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ERMS.integrated.markers, file.path(analysis_dir, "12_markers_SCT.csv"))

## Clean-up Seurat obj before saving
ERMS.integrated$nFeature_ADT <- NULL
ERMS.integrated$nCount_ADT <- NULL

## Save dataset
saveRDS(ERMS.integrated, file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ERMS_20230713.rds"))


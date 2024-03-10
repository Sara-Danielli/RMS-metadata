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

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/P3F1analysis')
if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## read file
P3F1.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_scores.rds"))
DefaultAssay(P3F1.integrated) <- "integrated"


## Order factors
P3F1.integrated$origin <- factor(x = P3F1.integrated$origin, levels = c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al."))
P3F1.integrated$model <- factor(x = P3F1.integrated$model, levels = c("Patient", "O-PDX", "Primary culture", "Cell line"))
P3F1.integrated$subtype <- factor(x = P3F1.integrated$subtype, levels = c("eRMS", "aRMS"))
P3F1.integrated$name <- factor(x = P3F1.integrated$name, levels = c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
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
P3F1.integrated <- CellCycleScoring(P3F1.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
P3F1.integrated <- ScaleData(P3F1.integrated)

## Define high-cycling vs low-cycling cells
Cycling <- WhichCells(P3F1.integrated, expression = S.Score > 0 | G2M.Score > 0)
nonCycling <- WhichCells(P3F1.integrated, expression = S.Score <= 0 & G2M.Score <= 0)

# Add cycling info to metadata column called 'Cycling_prop'
P3F1.integrated$Cycling_prop <- ifelse(colnames(P3F1.integrated) %in% Cycling, "Cycling", "Non-cycling")


## Run PCA
P3F1.integrated <- RunPCA(P3F1.integrated, reduction.name = "pca_rpca", reduction.key = "pca_rpca")

## Elbow plot to define number of PCs to use for downstream analysis    
ElbowPlot(P3F1.integrated, reduction = 'pca_rpca')
ggsave(file.path(analysis_dir, "1_Elbowplot_rpca.pdf"), width=6, height=4, dpi=300)

## Run UMAP and find clusters
P3F1.integrated <- FindNeighbors(P3F1.integrated, reduction = "pca_rpca", dims = 1:12)

# Define optimal number of clusters with clustree
PDX_clusters <- FindClusters(P3F1.integrated, resolution = seq(0, 0.5, 0.1))
clustree(PDX_clusters)
ggsave(file.path(analysis_dir,"2_clustre.pdf"), width=6, height=8, dpi=300)

P3F1.integrated <- FindClusters(P3F1.integrated, resolution = 0.3)

P3F1.integrated <- RunUMAP(P3F1.integrated, dims = 1:12, reduction = "pca_rpca",  reduction.key = "umap_rpca", reduction.name = "umap_rpca")


## change names of clusters
levels(P3F1.integrated)
new.cluster.ids <- c('DNA replication', 'Ground', 'Proliferative',  'Progenitor', 'Apoptosis', 'Differentiated', 'Neuronal')
names(new.cluster.ids) <- levels(P3F1.integrated)
P3F1.integrated <- RenameIdents(P3F1.integrated, new.cluster.ids)
P3F1.integrated[["cluster_names"]] <- Idents(object = P3F1.integrated)
Idents(P3F1.integrated) = "cluster_names"
P3F1.integrated$cluster_names <- factor(x = P3F1.integrated$cluster_names, 
                                        levels = c('Progenitor', 'Proliferative','DNA replication', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis'))

levels(P3F1.integrated)
new.cluster.ids.aggregate <- c('Progenitor', 'Proliferative', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')
names(new.cluster.ids.aggregate) <- levels(P3F1.integrated)
P3F1.integrated <- RenameIdents(P3F1.integrated, new.cluster.ids.aggregate)
P3F1.integrated[["cluster_names_aggregate"]] <- Idents(object = P3F1.integrated)
levels(P3F1.integrated) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')


Idents(P3F1.integrated) = P3F1.integrated$cluster_names

## Find markers
P3F1.integrated.markers <- FindAllMarkers(P3F1.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(P3F1.integrated.markers, file.path(analysis_dir, "12_markers_SCT.csv"))

## Clean-up Seurat obj before saving
P3F1.integrated$nFeature_ADT <- NULL
P3F1.integrated$nCount_ADT <- NULL
P3F1.integrated$GOBP_Neurogenesis1 <- NULL
P3F1.integrated$GOBP_Neuron_development1 <- NULL
P3F1.integrated$Wei_Neuronal_72117_21 <- NULL
P3F1.integrated$Wei_Neuronal_Mast1181 <- NULL
P3F1.integrated$Wei_Neuronal_Mast951 <- NULL
P3F1.integrated$Wei_Neuronal_747111 <- NULL
P3F1.integrated$Wei_Neuronal_72117_11 <- NULL

## Save dataset
saveRDS(P3F1.integrated, file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_20230713.rds"))


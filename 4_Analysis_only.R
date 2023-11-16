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

base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
data_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data/write')
plot_dir <- file.path(base_dir, 'output/metadata/Pseudobulk')

#base_dir <- "/n/scratch3/users/s/sad167/metadata"
source(file.path(base_dir, "Functions.R"))


#analysis_dir <- file.path(base_dir, 'analysis')
#if (!dir.exists(analysis_dir)){dir.create(analysis_dir, recursive = T)}

## read file
PDX.integrated <- readRDS(file.path(base_dir, "Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))
DefaultAssay(PDX.integrated) <- "integrated"


## Order factors
PDX.integrated$origin <- factor(x = PDX.integrated$origin, levels = c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al."))
PDX.integrated$model <- factor(x = PDX.integrated$model, levels = c("Patient", "O-PDX", "Primary culture", "Cell line"))
PDX.integrated$subtype <- factor(x = PDX.integrated$subtype, levels = c("eRMS", "aRMS"))
PDX.integrated$name <- factor(x = PDX.integrated$name, levels = c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
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
PDX.integrated <- CellCycleScoring(PDX.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## Define high-cycling vs low-cycling cells
Cycling <- WhichCells(PDX.integrated, expression = S.Score > 0 | G2M.Score > 0)
nonCycling <- WhichCells(PDX.integrated, expression = S.Score <= 0 & G2M.Score <= 0)

# Add cycling info to metadata column called 'Cycling_prop'
PDX.integrated$Cycling_prop <- ifelse(colnames(PDX.integrated) %in% Cycling, "Cycling", "Non-cycling")


## Run PCA
PDX.integrated <- RunPCA(PDX.integrated, reduction.name = "pca_rpca", reduction.key = "pca_rpca")

# ## Elbow plot to define number of PCs to use for downstream analysis    
ElbowPlot(PDX.integrated, reduction = 'pca_rpca')
ggsave(file.path(analysis_dir, "1_Elbowplot_rpca.pdf"), width=6, height=4, dpi=300)

## Run UMAP and find clusters
PDX.integrated <- FindNeighbors(PDX.integrated, reduction = "pca_rpca", dims = 1:6)

# Define optimal number of clusters with clustree
PDX_clusters <- FindClusters(PDX.integrated, resolution = seq(0, 0.5, 0.1))
clustree(PDX_clusters)
ggsave(file.path(analysis_dir,"2_clustre.pdf"), width=6, height=8, dpi=300)

PDX.integrated <- FindClusters(PDX.integrated, resolution = 0.3)

PDX.integrated <- RunUMAP(PDX.integrated, dims = 1:6, reduction = "pca_rpca",  reduction.key = "umap_rpca", reduction.name = "umap_rpca")


## change names of clusters
levels(PDX.integrated)

new.cluster.ids <- c('0-Ground',  '1-TR-progenitor',  '2-DNA replication', '3-Proliferative', '4-Ground', '5-DNA replication', '6-TR-differentiated', 
                     '7-Apoptosis', '8-Progenitor', '9-Differentiated', '10-DNA replication', '11-Ground')
names(new.cluster.ids) <- levels(PDX.integrated)
PDX.integrated <- RenameIdents(PDX.integrated, new.cluster.ids)
PDX.integrated[["cluster_names"]] <- Idents(object = PDX.integrated)
levels(PDX.integrated) <- c('8-Progenitor','1-TR-progenitor', '3-Proliferative', '2-DNA replication',  '5-DNA replication', '10-DNA replication',
                            '4-Ground','0-Ground', '11-Ground','6-TR-differentiated',  '9-Differentiated', '7-Apoptosis')

levels(PDX.integrated)
new.cluster.ids.aggregate <- c('Progenitor','TR-progenitor', 'Proliferative', 'Proliferative',  'Proliferative', 'Proliferative',
                               'Ground','Ground', 'Ground','TR-differentiated',  'Differentiated', 'Apoptosis')
names(new.cluster.ids.aggregate) <- levels(PDX.integrated)
PDX.integrated <- RenameIdents(PDX.integrated, new.cluster.ids.aggregate)
PDX.integrated[["cluster_names_aggregate"]] <- Idents(object = PDX.integrated)

Idents(PDX.integrated) = PDX.integrated$cluster_names

## Find markers
PDX.integrated.markers <- FindAllMarkers(PDX.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(PDX.integrated.markers, file.path(analysis_dir, "12_markers_SCT.csv"))

## Clean-up Seurat obj before saving
PDX.integrated$nFeature_ADT <- NULL
PDX.integrated$nCount_ADT <- NULL
PDX.integrated$GOBP_Neurogenesis1 <- NULL
PDX.integrated$GOBP_Neuron_development1 <- NULL
PDX.integrated$Wei_Neuronal_72117_21 <- NULL
PDX.integrated$Wei_Neuronal_Mast1181 <- NULL
PDX.integrated$Wei_Neuronal_Mast951 <- NULL
PDX.integrated$Wei_Neuronal_747111 <- NULL
PDX.integrated$Wei_Neuronal_72117_11 <- NULL
PDX.integrated$DS.Difference_Wei <- NULL
PDX.integrated$DS.Difference_Danielli <- NULL
PDX.integrated$DS.Difference_Patel <- NULL

## Save dataset
saveRDS(PDX.integrated, file.path(analysis_dir, "Danielli_Patel_Langenau_20230710.rds"))


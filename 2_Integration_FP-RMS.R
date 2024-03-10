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
plan("multicore", workers = 256)
options(future.globals.maxSize = 1000000 * 1024^2)

##################################################################
#  (1) Integration with RPCA correction ALL SAMPLES:  #
##################################################################

PDX.combined <- readRDS("/mnt/Sara/write/Danielli_Patel_Langenau_20221220.rds")

DefaultAssay(PDX.combined) <- "RNA"

## Subset subtypes
aRMS.combined <- subset(PDX.combined, subset = subtype == 'aRMS')

## Split object by name (= sample)
PDX.list <- SplitObject(aRMS.combined , split.by = "name")
PDX.list

PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# Select features for downstream integration
features <- SelectIntegrationFeatures(object.list = PDX.list)
PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
})


## Find anchors
anchors <- FindIntegrationAnchors(object.list = PDX.list, anchor.features = features, reduction = "rpca")

PDX.integrated <- IntegrateData(anchorset = anchors)


#  Scoring for cell cycle (has to be done before scaling data):  #
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
aRMS.combined <- CellCycleScoring(aRMS.combined, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(PDX.integrated) <- "integrated"


###### Save integrated dataset ###########
saveRDS(PDX.integrated, file = "/mnt/Sara/write/Danielli_Patel_Langenau_RPCA_ARMS_20230209.rds")



##################################################################
#  Scoring ARMS cells for gene sets:  #
##################################################################

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


PDX.integrated <- readRDS("/mnt/Sara/write/Danielli_Patel_Langenau_RPCA_ARMS_20230209.rds")
PDX.integrated 
DefaultAssay(PDX.integrated) <- "integrated"

############################################
############# Add scores ###################
############################################

## Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PDX.integrated <- CellCycleScoring(PDX.integrated, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## Add module scores
Wei_FN <- read_excel("/mnt/Sara/lists/Wei_2022/FN.xlsx", col_names=FALSE)
Wei_FN_genes <- as.list(Wei_FN)
Wei_FP <- read_excel("/mnt/Sara/lists/Wei_2022/FP.xlsx", col_names=FALSE)
Wei_FP_genes <- as.list(Wei_FP)
Wei_FP <- read_excel("/mnt/Sara/lists/Wei_2022/FP.xlsx", col_names=FALSE)
Wei_FP_genes <- as.list(Wei_FP)
Wei_Neuronal_Mast118 <- read_excel("/mnt/Sara/lists/Wei_2022/Neuronal_Mast118.xlsx", col_names=FALSE)
Wei_Neuronal_Mast118_genes <- as.list(Wei_Neuronal_Mast118)
Wei_Neuronal_Mast95 <- read_excel("/mnt/Sara/lists/Wei_2022/Neuronal_Mast95.xlsx", col_names=FALSE)
Wei_Neuronal_Mast95_genes <- as.list(Wei_Neuronal_Mast95)
Wei_Neuronal_74711 <- read_excel("/mnt/Sara/lists/Wei_2022/Neuronal_74711.xlsx", col_names=FALSE)
Wei_Neuronal_74711_genes <- as.list(Wei_Neuronal_74711)
Wei_Neuronal_72117_1 <- read_excel("/mnt/Sara/lists/Wei_2022/Neuronal_72117_4.xlsx", col_names=FALSE)
Wei_Neuronal_72117_1_genes <- as.list(Wei_Neuronal_72117_1)
Wei_Neuronal_72117_2 <- read_excel("/mnt/Sara/lists/Wei_2022/Neuronal_72117_3_5.xlsx", col_names=FALSE)
Wei_Neuronal_72117_2_genes <- as.list(Wei_Neuronal_72117_2)
Wei_Neuronal_72117_2 <- read_excel("/mnt/Sara/lists/Wei_2022/Neuronal_72117_3_5.xlsx", col_names=FALSE)
Wei_Neuronal_72117_2_genes <- as.list(Wei_Neuronal_72117_2)
GOBP_Neurogenesis <- read_excel("/mnt/Sara/lists/Wei_2022/GOBP_NEUROGENESIS.v2022.1.Hs.xlsx", col_names=FALSE)
GOBP_Neurogenesis_genes <- as.list(GOBP_Neurogenesis)
GOBP_Neuron_development <- read_excel("/mnt/Sara/lists/Wei_2022/GOBP_NEURON_DEVELOPMENT.v2022.1.Hs.xlsx", col_names=FALSE)
GOBP_Neuron_development_genes <- as.list(GOBP_Neuron_development)



Common_differentiated <- read_excel("/mnt/Sara/lists/Common_differentiated.xlsx", col_names=FALSE)
Common_differentiated <- as.list(Common_differentiated)
Common_proliferative <- read_excel("/mnt/Sara/lists/Common_proliferative.xlsx", col_names=FALSE)
Common_proliferative <- as.list(Common_proliferative)
Common_Stemcell <- read_excel("/mnt/Sara/lists/Common_Stemcell.xlsx", col_names=FALSE)
Common_Stemcell <- as.list(Common_Stemcell)


PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_Neuronal_Mast118_genes, ctrl = 100, name = "Wei_Neuronal_Mast118")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_Neuronal_Mast95_genes, ctrl = 100, name = "Wei_Neuronal_Mast95")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_Neuronal_74711_genes, ctrl = 100, name = "Wei_Neuronal_74711")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_Neuronal_72117_1_genes, ctrl = 100, name = "Wei_Neuronal_72117_1")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = GOBP_Neurogenesis_genes, ctrl = 100, name = "GOBP_Neurogenesis")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = GOBP_Neuron_development_genes, ctrl = 100, name = "GOBP_Neuron_development")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_Neuronal_72117_2_genes, ctrl = 100, name = "Wei_Neuronal_72117_2")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_differentiated, ctrl = 100, name = "Common_differentiated")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_Stemcell, ctrl = 100, name = "Common_Stemcell")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_proliferative, ctrl = 100, name = "Common_proliferative")

## Calculate stemness score (= stem-like - differentiated)
PDX.integrated$DS.Difference_common <- PDX.integrated$Common_differentiated1 - PDX.integrated$Common_Stemcell1


# Scale data 
PDX.integrated <- ScaleData(PDX.integrated, verbose = TRUE)


###### Save integrated dataset ###########
saveRDS(PDX.integrated, file = "/mnt/Sara/write/Danielli_Patel_Langenau_RPCA_ARMS_scores_20230209.rds")


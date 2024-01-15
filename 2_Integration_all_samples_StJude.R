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
plan("multicore", workers = 128)
options(future.globals.maxSize = 1000000 * 1024^2)

##################################################################
#  (1) Integration with RPCA correction ALL SAMPLES:  #
##################################################################

PDX.combined <- readRDS("./Danielli_Patel_Langenau_20221220.rds")

DefaultAssay(PDX.combined) <- "RNA"

## Split object by name (= sample)
PDX.list <- SplitObject(PDX.combined , split.by = "name")
PDX.list

PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

saveRDS(PDX.list, file = "./Danielli_Patel_Langenau_20221220_norm.rds")

# Select features for downstream integration
features <- SelectIntegrationFeatures(object.list = PDX.list)
PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
})

saveRDS(PDX.list, file = "./Danielli_Patel_Langenau_20221220_scale.rds")

## Find anchors
anchors <- FindIntegrationAnchors(object.list = PDX.list, anchor.features = features, reduction = "rpca")

PDX.integrated <- IntegrateData(anchorset = anchors)


#  Scoring for cell cycle (has to be done before scaling data):  #
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PDX.combined <- CellCycleScoring(PDX.combined, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(PDX.integrated) <- "integrated"


###### Save integrated dataset ###########
saveRDS(PDX.integrated, file = "./Danielli_Patel_Langenau_RPCA_20221221.rds")



##################################################################
#  (2) Integration with RPCA correction aRMS:  #
##################################################################

PDX.combined <- readRDS("./Danielli_Patel_Langenau_20221220.rds")

## Subset subtypes
aRMS.combined <- subset(aRMS.combined, subset = subtype == 'aRMS')

DefaultAssay(aRMS.combined) <- "RNA"

## Split object by name (= sample)
PDX.list <- SplitObject(aRMS.combined , split.by = "name")
PDX.list

PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

saveRDS(PDX.list, file = "./Danielli_Patel_Langenau_aRMS_norm.rds")

# Select features for downstream integration
features <- SelectIntegrationFeatures(object.list = PDX.list)
PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
})

saveRDS(PDX.list, file = "./Danielli_Patel_Langenau_aRMS_scale.rds")

## Find anchors
anchors <- FindIntegrationAnchors(object.list = PDX.list, anchor.features = features, reduction = "rpca")

aRMS.integrated <- IntegrateData(anchorset = anchors)


#  Scoring for cell cycle (has to be done before scaling data):  #
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
aRMS.combined <- CellCycleScoring(aRMS.combined, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(aRMS.integrated) <- "integrated"


###### Save integrated dataset ###########
saveRDS(aRMS.integrated, file = "./Danielli_Patel_Langenau_RPCA_aRMS.rds")




##################################################################
#  (3) Integration with RPCA correction eRMS:  #
##################################################################

PDX.combined <- readRDS("./Danielli_Patel_Langenau_20221220.rds")

## Subset subtypes
eRMS.combined <- subset(eRMS.combined, subset = subtype == 'eRMS')

DefaultAssay(eRMS.combined) <- "RNA"

## Split object by name (= sample)
PDX.list <- SplitObject(eRMS.combined , split.by = "name")
PDX.list

PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

saveRDS(PDX.list, file = "./Danielli_Patel_Langenau_eRMS_norm.rds")

# Select features for downstream integration
features <- SelectIntegrationFeatures(object.list = PDX.list)
PDX.list <- lapply(X = PDX.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
})

saveRDS(PDX.list, file = "./Danielli_Patel_Langenau_eRMS_scale.rds")

## Find anchors
anchors <- FindIntegrationAnchors(object.list = PDX.list, anchor.features = features, reduction = "rpca")

eRMS.integrated <- IntegrateData(anchorset = anchors)


#  Scoring for cell cycle (has to be done before scaling data):  #
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
eRMS.combined <- CellCycleScoring(eRMS.combined, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(eRMS.integrated) <- "integrated"


###### Save integrated dataset ###########
saveRDS(eRMS.integrated, file = "./Danielli_Patel_Langenau_RPCA_eRMS.rds")




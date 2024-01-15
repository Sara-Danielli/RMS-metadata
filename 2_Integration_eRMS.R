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

base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
genelist_dir <- file.path(base_dir, 'lists')

##################################################################
#  (1) Integration with RPCA correction ERMS:  #
##################################################################

PDX.combined <- readRDS("/mnt/Sara/write/Danielli_Patel_Langenau_20221220.rds")

DefaultAssay(PDX.combined) <- "RNA"

## Subset subtypes
eRMS.combined <- subset(PDX.combined, subset = subtype == 'eRMS')

## Split object by name (= sample)
PDX.list <- SplitObject(eRMS.combined , split.by = "name")
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
eRMS.combined <- CellCycleScoring(eRMS.combined, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(PDX.integrated) <- "integrated"


###### Save integrated dataset ###########
saveRDS(PDX.integrated, file = "/mnt/Sara/write/Danielli_Patel_Langenau_RPCA_ERMS_20230215.rds")




##################################################################
#  Scoring ERMS cells for gene sets:  #
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


PDX.integrated <- readRDS("/mnt/Sara/write/Danielli_Patel_Langenau_RPCA_ERMS_20230215.rds")
PDX.integrated 
DefaultAssay(PDX.integrated) <- "integrated"

############################################
############# Add scores ###################
############################################

## Add module scores
Common_differentiated <- read_excel("/mnt/Sara/lists/Common_differentiated.xlsx", col_names=FALSE)
Common_differentiated <- as.list(Common_differentiated)
Common_proliferative <- read_excel(file.path(genelist_dir, "Common_proliferative.xlsx"), col_names=FALSE)
Common_proliferative <- as.list(Common_proliferative)
Common_Stemcell <- read_excel("/mnt/Sara/lists/Common_Stemcell.xlsx", col_names=FALSE)
Common_Stemcell <- as.list(Common_Stemcell)


PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_differentiated, ctrl = 100, name = "Common_differentiated")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_Stemcell, ctrl = 100, name = "Common_Stemcell")
PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_proliferative, ctrl = 100, name = "Common_proliferative")

## Calculate stemness score (= stem-like - differentiated)
PDX.integrated$DS.Difference_common <- PDX.integrated$Common_differentiated1 - PDX.integrated$Common_Stemcell1

## Cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
PDX.integrated <- CellCycleScoring(PDX.integrated, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Scale data 
PDX.integrated <- ScaleData(PDX.integrated, verbose = TRUE)


###### Save integrated dataset ###########
saveRDS(PDX.integrated, file.path(base_dir,"write/Danielli_Patel_Langenau_RPCA_ERMS_scores_20230215.rds"))
PDX.integrated <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ERMS_scores_20230215.rds"))




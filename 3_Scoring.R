rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(data.table)
library(SeuratObject)

base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data/metadata')
#base_dir <- file.path("/n/scratch3/users/s/sad167/HOPE")
genelist_dir <- file.path(base_dir, 'lists')

write_dir <- file.path(base_dir, 'write')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}

# read integrated object
#PDX.integrated <- readRDS(file.path(base_dir, "Danielli_Patel_Langenau_RPCA_20221221.rds"))
PDX.integrated <- readRDS(file.path(base_dir, "Danielli_Patel_Langenau_20230710.rds"))

PDX.integrated 
DefaultAssay(PDX.integrated) <- "integrated"


# Load gene signatures

        ## Cell cycle scoring
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
        PDX.integrated <- CellCycleScoring(PDX.integrated, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

        ## Publication-specific scores
        Wei_mesenchymal <- read_excel(file.path(genelist_dir, "Wei_2022/Mesenchymal.xlsx", col_names=FALSE))
        Wei_mesenchymal_genes <- as.list(Wei_mesenchymal)
        Wei_proliferative <- read_excel(file.path(genelist_dir, "Wei_2022/Proliferative.xlsx", col_names=FALSE))
        Wei_proliferative_genes <- as.list(Wei_proliferative)
        Wei_muscle <- read_excel(file.path(genelist_dir, "Wei_2022/Muscle.xlsx", col_names=FALSE))
        Wei_muscle_genes <- as.list(Wei_muscle)
        Wei_FN <- read_excel(file.path(genelist_dir, "Wei_2022/FN.xlsx", col_names=FALSE))
        Wei_FN_genes <- as.list(Wei_FN)
        Wei_FP <- read_excel(file.path(genelist_dir, "Wei_2022/FP.xlsx", col_names=FALSE))
        Wei_FP_genes <- as.list(Wei_FP)
        Wei_FP <- read_excel(file.path(genelist_dir, "Wei_2022/FP.xlsx", col_names=FALSE))
        Wei_FP_genes <- as.list(Wei_FP)

        Danielli_MuSC <- read_excel(file.path(genelist_dir, "Danielli_2022/MuSClike.xlsx", col_names=FALSE))
        Danielli_MuSC_genes <- as.list(Danielli_MuSC)
        Danielli_cycling <- read_excel(file.path(genelist_dir, "Danielli_2022/Cycling.xlsx", col_names=FALSE))
        Danielli_cycling_genes <- as.list(Danielli_cycling)
        Danielli_differentiated <- read_excel(file.path(genelist_dir, "Danielli_2022/Differentiated.xlsx", col_names=FALSE))
        Danielli_differentiated_genes <- as.list(Danielli_differentiated)
        
        Patel_mesoderm <- read_excel(file.path(genelist_dir, "Patel_2022/Mesoderm.xlsx", col_names=FALSE))
        Patel_mesoderm_genes <- as.list(Patel_mesoderm)
        Patel_myoblast <- read_excel(file.path(genelist_dir, "Patel_2022/Myoblast.xlsx", col_names=FALSE))
        Patel_myoblast_genes <- as.list(Patel_myoblast)
        Patel_myocyte <- read_excel(file.path(genelist_dir, "Patel_2022/Myocyte.xlsx", col_names=FALSE))
        Patel_myocyte_genes <- as.list(Patel_myocyte)
        
        ## RMS atlas scores (identified in this new analysis)
        Common_differentiated <- read_excel(file.path(genelist_dir, "Common_differentiated.xlsx", col_names=FALSE))
        Common_differentiated <- as.list(Common_differentiated)
        Common_proliferative <- read_excel(file.path(genelist_dir,"Common_proliferative.xlsx", col_names=FALSE))
        Common_proliferative <- as.list(Common_proliferative)
        Common_Stemcell <- read_excel(file.path(genelist_dir, "Common_Stemcell.xlsx", col_names=FALSE))
        Common_Stemcell <- as.list(Common_Stemcell)
        
        
# Add scores
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_mesenchymal_genes, ctrl = 100, name = "Wei_mesenchymal")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_proliferative_genes, ctrl = 100, name = "Wei_proliferative")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_muscle_genes, ctrl = 100, name = "Wei_muscle")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_FP_genes, ctrl = 100, name = "Wei_FP")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Wei_FN_genes, ctrl = 100, name = "Wei_FN")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Danielli_MuSC_genes, ctrl = 100, name = "Danielli_MuSC")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Danielli_cycling_genes, ctrl = 100, name = "Danielli_cycling")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Danielli_differentiated_genes, ctrl = 100, name = "Danielli_differentiated")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Patel_mesoderm_genes, ctrl = 100, name = "Patel_mesoderm")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Patel_myoblast_genes, ctrl = 100, name = "Patel_myoblast")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Patel_myocyte_genes, ctrl = 100, name = "Patel_myocyte")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_differentiated, ctrl = 100, name = "Common_differentiated")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_Stemcell, ctrl = 100, name = "Common_Stemcell")
        PDX.integrated <- AddModuleScore(object = PDX.integrated, assay = 'RNA', features = Common_proliferative, ctrl = 100, name = "Common_proliferative")

# Calculate muscle lineage score (DS.difference = differentiated - progenitor)
        PDX.integrated$DS.Difference_Wei <- PDX.integrated$Wei_muscle1 - PDX.integrated$Wei_mesenchymal1
        PDX.integrated$DS.Difference_Danielli <- PDX.integrated$Danielli_differentiated1 - PDX.integrated$Danielli_MuSC1
        PDX.integrated$DS.Difference_Patel <- PDX.integrated$Patel_myocyte1 - PDX.integrated$Patel_mesoderm1
        PDX.integrated$DS.Difference_common <- PDX.integrated$Common_differentiated1 - PDX.integrated$Common_Stemcell1
        
# Add Cell cycle score
        s.genes <- cc.genes$s.genes
        g2m.genes <- cc.genes$g2m.genes
        PDX.integrated <- CellCycleScoring(PDX.integrated, assay = 'RNA', s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        
# Scale data 
        PDX.integrated <- ScaleData(PDX.integrated, verbose = TRUE)


# Save scored dataset
#saveRDS(PDX.integrated, file = "/mnt/Sara/write/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds")
saveRDS(PDX.integrated, file.path(write_dir, "Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))
        



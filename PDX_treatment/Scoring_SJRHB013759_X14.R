rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(data.table)
library(SingleCellExperiment)
library(paletteer)
library(SeuratObject)
library(DESeq2)
library(scater)
library(tidyverse)
library(cowplot)
library(Matrix.utils)
library(edgeR)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(apeglm)
library(DESeq2)
library(UCell)
library(ggpubr)
plan("multiprocess", workers = 256)
options(future.globals.maxSize = 256000 * 1024^2)



##################################################################
#  (0) Loading Patel FFPE RNAseq datasets:  #
##################################################################

### ! Excel matrix contained some numbers that were not integers --> approximated

# load TPM data
dataset <- read.delim("/mnt/Sara/RMS/Patel_2022/FFPE_RNAseq/aRMS/log2_data_SJRHB013759_X14.txt")
    # set gene names as row header
    dataset3 <- dataset[, -1]
    
    rownames(dataset3) <- make.names(dataset$geneSymbol, unique = TRUE)
      
# load metadata 
FFPE_meta <- read_excel("/mnt/Sara/RMS/Patel_2022/FFPE_RNAseq/aRMS/metadata_log2_data_SJRHB013759_X14.xlsx", col_names=TRUE)
      # set first column as row header
    FFPE_meta2 <- FFPE_meta %>% remove_rownames %>% column_to_rownames(var="id")

# Create Seurat object
FFPE_seurat <- CreateSeuratObject(counts = dataset3, project = "FFPE")
FFPE_seurat <- AddMetaData(FFPE_seurat, metadata = FFPE_meta2)

##################################################################
#  (1) Score datasets with gene signatures:  #
##################################################################

# Load signatures identified in the integrated dataset
PAX3FOXO1_markers <- read_excel("/mnt/Sara/lists/FFPE_scoring/PAX3FOXO1_markers.xlsx", col_names=FALSE)

  # Set column name
  colnames(PAX3FOXO1_markers) <- PAX3FOXO1_markers[1,]
  
  # remove column name
  PAX3FOXO1_markers <- PAX3FOXO1_markers[-1,]
  
  neuronal <- PAX3FOXO1_markers$neuronal[1:288]
  progenitor <- PAX3FOXO1_markers$progenitor[1:185]
  differentiated <- PAX3FOXO1_markers$differentiated[1:257]
  proliferative <- PAX3FOXO1_markers$proliferative[1:158]
  DNA_replication <- PAX3FOXO1_markers$DNA_replication[1:118]
  

# Add scores using "UCell" fct
    markers <- list()
    markers$differentiated <- differentiated
    markers$proliferative <- proliferative
    markers$progenitor <- progenitor
    markers$neuronal <- neuronal
    markers$DNA_replication <- DNA_replication
    
    markers$Common_differentiated <- read_excel("/mnt/Sara/lists/Common_differentiated.xlsx", col_names=FALSE)
    markers$Common_proliferative <- read_excel("/mnt/Sara/lists/Common_proliferative.xlsx", col_names=FALSE)
    markers$Common_Stemcell <- read_excel("/mnt/Sara/lists/Common_Stemcell.xlsx", col_names=FALSE)
    
    
    set.seed(123)
    FFPE_seurat <- AddModuleScore_UCell(FFPE_seurat, features = markers)
    signature.names <- paste0(names(markers), "_UCell")
    
    
# Add scores using "Addmodulescore" fct
    neuronal <- list(neuronal)
    progenitor <- list(progenitor)
    differentiated <- list(differentiated)
    proliferative <- list(proliferative)
    DNA_replication <- list(DNA_replication)
    
    FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = differentiated, ctrl = 100, name = "differentiated")
    FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = proliferative, ctrl = 100, name = "proliferative")
    FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = progenitor, ctrl = 100, name = "progenitor")
    FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = neuronal, ctrl = 100, name = "neuronal")
    FFPE_seurat <- AddModuleScore(object = FFPE_seurat, assay = 'RNA', features = DNA_replication, ctrl = 100, name = "DNA_replication")
    


# Scale data
ScaleData(FFPE_seurat)



## Export data tables
md <- FFPE_seurat@meta.data %>% as.data.table
table_PDX <- md[, .N, by = c("name", "neuronal_UCell", "DNA_replication_UCell", "progenitor_UCell", "differentiated_UCell", "proliferative_UCell", 
                             "Common_differentiated_UCell", "Common_proliferative_UCell", "Common_Stemcell_UCell")] 
write.csv(table_PDX, "/mnt/Sara/output/metadata/FFPE_Patel/SJRHB013759_X14/Signatures_UCell.csv")

table_PDX <- md[, .N, by = c("name", "differentiated1", "proliferative1", "progenitor1", "neuronal1", "DNA_replication1")] 
write.csv(table_PDX, "/mnt/Sara/output/metadata/FFPE_Patel/SJRHB013759_X14/Signatures_Addmodulescore.csv")

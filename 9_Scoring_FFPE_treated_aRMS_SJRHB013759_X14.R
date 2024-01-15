rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(data.table)
library(paletteer)
library(SeuratObject)
library(tidyverse)
library(cowplot)
library(ggpubr)

base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
data_dir <- file.path(base_dir, 'RMS/Patel_2022/FFPE_RNAseq/aRMS')
plot_dir <- file.path(base_dir, 'output/metadata/Pseudobulk/PDX')

if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}




##################################################################
#  (0) Loading Patel FFPE RNAseq datasets:  #
##################################################################

### ! Excel matrix contained some numbers that were not integers --> approximated

# load TPM data
dataset <- read.delim(file.path(data_dir, "log2_data_SJRHB013759_X14.txt"))
    # set gene names as row header
    dataset3 <- dataset[, -1]
    
    rownames(dataset3) <- make.names(dataset$geneSymbol, unique = TRUE)
      
# load metadata 
FFPE_meta <- read_excel(file.path(data_dir, "metadata_log2_data_SJRHB013759_X14.xlsx"), col_names=TRUE)
      # set first column as row header
    FFPE_meta2 <- FFPE_meta %>% remove_rownames %>% column_to_rownames(var="id")

# Create Seurat object
FFPE_seurat <- CreateSeuratObject(counts = dataset3, project = "FFPE")
FFPE_seurat <- AddMetaData(FFPE_seurat, metadata = FFPE_meta2)

##################################################################
#  (1) Score datasets with gene signatures:  #
##################################################################

# Load signatures identified in the integrated dataset
PAX3FOXO1_markers <- read_excel("/Volumes/Sara_PhD/scRNAseq_data/lists/FFPE_scoring/PAX3FOXO1_markers.xlsx", col_names=FALSE)

  # Set column name
  colnames(PAX3FOXO1_markers) <- PAX3FOXO1_markers[1,]
  
  # remove column name
  PAX3FOXO1_markers <- PAX3FOXO1_markers[-1,]
  
  neuronal <- PAX3FOXO1_markers$neuronal[1:288]
  progenitor <- PAX3FOXO1_markers$progenitor[1:185]
  differentiated <- PAX3FOXO1_markers$differentiated[1:257]
  proliferative <- PAX3FOXO1_markers$proliferative[1:259]


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

table_PDX <- md[, .N, by = c("name", "differentiated1", "proliferative1", "progenitor1", "neuronal1")] 
write.csv(table_PDX, file.path(plot_dir, "Table_Signatures_Addmodulescore.csv"))

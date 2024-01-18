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
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
data_dir <- file.path(base_dir, 'RMS/Patel_2022/FFPE_RNAseq')
plot_dir <- file.path(base_dir, 'output/metadata/Pseudobulk')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}

genelist_dir <- file.path(base_dir, 'list_final')


##################################################################
#  (0) Loading Patel FFPE RNAseq datasets:  #
##################################################################

## ! Excel matrix contained some numbers that were not integers --> approximated

# load TPM data
dataset <- read.delim(file.path(data_dir, "RMS13_FFPE_RNA-seq_tpm_data.txt"))

rownames(dataset) <- make.names(dataset$geneSymbol, unique = TRUE)

# set gene names as row header
dataset2 <- dataset[, -c(1:4)]

# log-normalize
dataset_log <-log1p(dataset2+1)

# subset to ARMS
#FPRMS <- dataset_log[ ,c('X19A', 'X19B', 'X21A', 'X21B')]
FPRMS <- dataset_log


# load metadata  -----------------------------------
FFPE_meta <- read_excel(file.path(data_dir,"metadata.xlsx"), col_names=TRUE)
  # set first column as row header
  FFPE_meta2 <- FFPE_meta %>% remove_rownames %>% column_to_rownames(var="patient_id")
# subset to FPRMS samples
  #metadata <- FFPE_meta2[c('X19A', 'X19B', 'X21A', 'X21B'), ]
  metadata <- FFPE_meta2
  


# subset to synaptic genes -----------------------------------
  # load gene list integrated dataset
  markers <- read_xlsx(file.path(genelist_dir, 'all_integrated_cluster_markers.xlsx'))
  
  # order by Annotation and fold change
  markers <- markers %>% arrange(Annotation, 'Average log2 fold change')
  markers <- as.data.frame(markers)
  
  # convert into list
  marker_list <- split(markers[, -c(2:7)], f = markers$Annotation)
  
 #  # load gene lists P3F1
 #  markers <- read.csv(file.path(genelist_dir, 'ERMS_cluster_markers.csv'))
 # 
 #  # rename DNA replication as Proliferative
 #  markers <- markers %>%
 #    mutate_all(~ ifelse(. == 'DNA replication', 'Proliferative', .)) %>%
 #    mutate_all(~ ifelse(. == '4-Progenitor', 'Progenitor', .)) %>%
 #    mutate_all(~ ifelse(. == '6-Progenitor', 'Progenitor', .)) %>%
 #    mutate_all(~ ifelse(. == '2-Ground', 'Ground', .)) 
 #  
 #  # order by Annotation and fold change
 #  markers <- markers %>% arrange(cluster, desc(avg_log2FC))
 #  markers <- as.data.frame(markers)
 #  
 #  # convert into list
 #  marker_list <- split(markers[, -c(1:7)], f = markers$cluster)
 # # marker_list <- lapply(marker_list, head, 50)


# Subset FFPE dataset to genes of interest -----------------------------------
  # create df with genes of interest
  genes <- marker_list
  genes_df <- data.frame(
    Values = unlist(marker_list),
    Group = rep(names(marker_list), sapply(marker_list, length))
  )
  rownames(genes_df) <- NULL
  
  seurat_df <- FPRMS[genes_df$Values, ]

  
# define annotations -----------------------------------
  # order dataframe based on clusters and then shuffle cells inside (random shuffling of Ds.difference)
  order <- data.frame(metadata %>%
                        arrange(name))
  ordered_rows <- rownames(order)
  
  # scale dataset
  seurat_df <- seurat_df - rowMeans(seurat_df)
  seurat_df <- as.data.frame(seurat_df)
  
  # order based on  score
  seurat_df <- seurat_df[, ordered_rows]
  
  # define annotations (important: order based on new order of cells)
  metadata2 <- metadata[colnames(seurat_df), ]
  metadata <- metadata2
  annotation_info <- metadata[, c("subtype", "treatment", "name")]
  
  
  # define colors  -----------------------------------
  # Color Louvain clusters
  col_Cluster_assignment <-  c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF' , '#D8D8E0FF')
  names(col_Cluster_assignment) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')
  
  col_subtype <- c('#D3A2C2FF', '#95CECFFF')
  names(col_subtype) <- c('eRMS', 'aRMS')
  
  col_column_ha <- list(subtype = col_subtype,
                        Cluster_assignment = col_Cluster_assignment)
  
  
  # define column annotation
  column_ha = HeatmapAnnotation(df = annotation_info
                                #col = col_column_ha
                                )
  
  # define row annotation
  #row_ha = rowAnnotation(Group = genes_df$Group)
  
  # genes to mark
  # elements_to_find <- c('FN1', 'CAV1',  'MYOD1', 'MYOG', 'MYH3', 'MKI67', 'CENPF', 'TTN', 'MYL4', 'DCX', 'SYT1', 'L1CAM', 'STMN4')
  # rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
  # rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]
  # 
  # right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
  #                                          labels = rows_genes_to_mark_name))
  
  
  # plot heatmap 
  Cairo::CairoPDF(file.path(plot_dir, "1_heatmap_integrated_cluster_markers.pdf"), width=10, height=10)
  ht <- Heatmap(seurat_df,
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                #bottom_annotation = bottom_ha,
                border = TRUE,
                top_annotation = column_ha,
                row_split = factor(c(genes_df$Group),levels = c('Progenitor', 'TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated',
                                                                 'Differentiated', 'Apoptosis')),
                #right_annotation = right_ha,
                use_raster = TRUE,
                raster_quality = 10,
                cluster_rows = FALSE
  )
  draw(ht)
  dev.off()
  
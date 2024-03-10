# Load packages -----------------------------------
rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(data.table)
library(ComplexHeatmap)
#library(magick)
#library(RColorBrewer)
library(circlize)
#library(ggpubr)
#library(openxlsx)

# Organize environment  -----------------------------------
#base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'
base_dir <- '/n/scratch/users/s/sad167/RMS'

genelist_dir <- file.path(base_dir, 'data/list_final')

analysis_dir <- file.path(base_dir, 'analysis/heatmap')
#analysis_dir <- file.path(base_dir, 'heatmap')

plot_dir <- file.path(analysis_dir, 'plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
write_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}



# -----------------------------------------------------------------------
# (4) Heatmap top marker genes FN-RMS clusters
# -----------------------------------------------------------------------
# load Seurat object -----------------------------------
ERMS <- readRDS(file.path(base_dir, "data/Danielli_Patel_Langenau_RPCA_ERMS_20230713.rds"))
# rename ERMS (problem with current IDs)
Idents(ERMS) = "cluster_names"
new.cluster.ids.aggregate <- c('Progenitor', 'Progenitor', 'Proliferative',  'Differentiated', 'IFN',  'Proliferative', 'Ground', 'Ground')
names(new.cluster.ids.aggregate) <- levels(ERMS)
ERMS <- RenameIdents(ERMS, new.cluster.ids.aggregate)
ERMS[["cluster_names_aggregate"]] <- Idents(object = ERMS)

#reorder
ERMS$cluster_names_aggregate <- factor(x = ERMS$cluster_names_aggregate , 
                                       levels = c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'IFN'))
# rename subtype
metadata <- ERMS@meta.data

# rename DS.Difference_common score into lineage score
colnames(ERMS@meta.data)[32] <- "cluster_names_aggregate"


# save metadata as df
metadata <- data.frame(ERMS@meta.data)

# aggregate cells from same cluster -----------------------------------
Idents(ERMS) <- 'cluster_names_aggregate'

ERMS <- AggregateExpression(ERMS, return.seurat = TRUE)

# load gene lists -----------------------------------
markers <- read.csv(file.path(genelist_dir, 'ERMS_cluster_markers.csv'))

# rename DNA replication as Proliferative
markers <- markers %>%
  mutate_all(~ ifelse(. == 'DNA replication', 'Proliferative', .)) %>%
  mutate_all(~ ifelse(. == '4-Progenitor', 'Progenitor', .)) %>%
  mutate_all(~ ifelse(. == '6-Progenitor', 'Progenitor', .))  %>%
  mutate_all(~ ifelse(. == '2-Ground', 'Ground', .))

# order by Annotation and fold change
markers <- markers %>% arrange(cluster, desc(avg_log2FC))
markers <- as.data.frame(markers)

# convert into list
marker_list <- split(markers[, -c(1:7)], f = markers$cluster)

# Remove duplicates within each vector
marker_list <- lapply(marker_list, unique)

# select top 100 genes per cluster
#marker_list <- lapply(marker_list,head,100)




# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- marker_list
genes_df <- data.frame(
  Values = unlist(marker_list),
  Group = rep(names(marker_list), sapply(marker_list, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- ERMS@assays$RNA$data[genes_df$Values, ]


# define annotations -----------------------------------

# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)


# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)

# genes to mark
elements_to_find <- c('FN1', 'CAV1', 'PAX7', 'MEOX2', 'CD44', 
                      'MYOD1', 'MYOG', 'MYH3', 
                      'MKI67', 'CENPF', 
                      'TTN', #'MYL4', 
                      'SYP', 'L1CAM', 'CHGA' # #'DCX', 'STMN4'
)

rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1, 0, 1), c("#2E5A87FF", "#FCFDFEFF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_ERMS_heatmap_aggregate_all.pdf"), width=5, height=6)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              #bottom_annotation = bottom_ha,
              border = TRUE,
              row_split = factor(c(genes_df$Group), 
                                 levels = c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'IFN')),
              #top_annotation = column_ha,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun,
              row_names_rot = 30
)
draw(ht)
dev.off()


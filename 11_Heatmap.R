# Load packages -----------------------------------
rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(data.table)
library(ComplexHeatmap)
library(magick)
library(RColorBrewer)
library(circlize)
library(ggpubr)
#library(openxlsx)

# Organize environment  -----------------------------------
base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'

genelist_dir <- file.path(base_dir, 'list_final')

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/heatmap')

plot_dir <- file.path(analysis_dir, 'plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
write_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}


# -----------------------------------------------------------------------
# (1) Heatmap combined RMS dataset scored for muscle lineage score
# -----------------------------------------------------------------------

# load Seurat object -----------------------------------
P3F1 <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_20230202_scoring100.rds"))

# rename subtype
  metadata <- PDX.combined@meta.data
  # Use dplyr's mutate to substitute specific cells with new names
  metadata <- metadata %>%
    mutate(subtype = case_when(
      subtype == "eRMS" ~ "FN-RMS",
      subtype == "aRMS" ~ "FP-RMS"
    ))
  
  # Add new metadata to Seurat objetct
  PDX.combined <- AddMetaData(PDX.combined, metadata)


# rename DS.Difference_common score into lineage score
colnames(PDX.combined@meta.data)[39] <- "Muscle_lineage_score"
colnames(PDX.combined@meta.data)[43] <- "Cluster_assignment"

# save metadata as df
metadata <- data.frame(PDX.combined@meta.data)

# load gene lists -----------------------------------
markers <- read_xlsx(file.path(genelist_dir, 'all_integrated_cluster_markers.xlsx'))

# order by Annotation and fold change
markers <- markers %>% arrange(Annotation, 'Average log2 fold change')
markers <- as.data.frame(markers)

# convert into list
marker_list <- split(markers[, -c(2:7)], f = markers$Annotation)

# select top 10 genes per cluster
marker_list <- lapply(marker_list,head,50)


# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- marker_list
genes_df <- data.frame(
  Values = unlist(marker_list),
  Group = rep(names(marker_list), sapply(marker_list, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- PDX.combined@assays$RNA@data[genes_df$Values, ]


# define annotations -----------------------------------
order <- data.frame(metadata %>% 
                      arrange(Muscle_lineage_score))
ordered_rows <- rownames(order)
ordered_scores <- order$Muscle_lineage_score


# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)

# order based on  score
seurat_df <- seurat_df[, ordered_rows]

# define annotations (important: order based on new order of cells)
metadata2 <- metadata[colnames(seurat_df), ]
metadata <- metadata2
annotation_info <- metadata[, c(
  #"name", 
  "subtype", "model", "Cluster_assignment", "Cycling_prop", "Muscle_lineage_score")]


# define colors  -----------------------------------
# Color model
col_model <- c("#009B9EFF","#A7D3D4FF",  "#E4C1D9FF","#C75DABFF")
names(col_model) <- c("Patient", "O-PDX", "Primary culture", "Cell line")

# Color origin
col_origin <- paletteer::paletteer_d("ggthemes::excel_Slice", n=4)
names(col_origin) <- c("Wei et al.", "Patel et al.", "Danielli et al.", "Weng et al.")

# Color name
col_aRMS <- paletteer::paletteer_c("ggthemes::Blue-Teal", n = 27)
col_eRMS <- paletteer::paletteer_c("ggthemes::Purple", n = 45)
col_name <- c(col_aRMS, col_eRMS)
names(col_name) <- c('20082',  'aRMS-1',  'aRMS-2', 'aRMS-3', 
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
                     'SJRHB030680_X1', 'SJRHB049189_D1', 'SJRHB049189_X1')

# Color subtype
col_subtype <- c('#D3A2C2FF', '#95CECFFF')
names(col_subtype) <- c('FN-RMS', 'FP-RMS')

# Color Louvain clusters
col_Cluster_assignment <- c("#276419FF","#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#DE77AEFF", "#8E0152FF", '#D8D8E0FF')
names(col_Cluster_assignment) <- c('Progenitor','TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated', 'Differentiated', 'Apoptosis')

# Color cycling properties
col_cycling <- c('black', 'grey')
names(col_cycling) <- c("Cycling", "Non-cycling")

col_scores <- colorRamp2(c(-1, 0, 1.5), c("#4d9221", "#f7f7f7", "#c51b7d"))


col_column_ha <- list(
  #name = col_name,
                      subtype = col_subtype,
                      model = col_model,
                      Cluster_assignment = col_Cluster_assignment,
                      Cycling_prop = col_cycling,
                      Muscle_lineage_score = col_scores
                      )


# define column annotation
column_ha = HeatmapAnnotation(df = annotation_info,
                              col = col_column_ha)

# define bottom annotation
bottom_ha = HeatmapAnnotation(scores = anno_barplot(ordered_scores), height = unit(1.5, "cm"))

# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)

# genes to mark
elements_to_find <- c('CD44', 'FN1', 'EGFR', 'MYOD1', 'MYOG', 'MYH3', 'MKI67', 'HELLS', 'TTN', 'MYL4')
rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1, -0.1, 1), c("#2E5A87FF", "#CCCCD2FF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_heatmap_50.pdf"), width=15, height=6)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              bottom_annotation = bottom_ha,
              border = TRUE,
              column_split = factor(c(annotation_info$subtype), levels = c('FN-RMS', 'FP-RMS')), 
              row_split = factor(c(genes_df$Group), 
                                 levels = c('Progenitor', 'TR-progenitor', 'Proliferative', 'Ground', 'TR-differentiated',
                                            'Differentiated', 'Apoptosis')),
              top_annotation = column_ha,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun
              )
draw(ht)
dev.off()




# -----------------------------------------------------------------------
# (2) Heatmap top marker genes FP-RMS (PAX3::FOXO1) clusters
# -----------------------------------------------------------------------

# load Seurat object -----------------------------------
P3F1 <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P3F1_20230713.rds"))
  # rename P3F1 (problem with current IDs)
  Idents(P3F1) = "cluster_names"
  new.cluster.ids.aggregate <- c('Progenitor', 'Proliferative', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')
  names(new.cluster.ids.aggregate) <- levels(P3F1)
  P3F1<- RenameIdents(P3F1, new.cluster.ids.aggregate)
  P3F1[["cluster_names_aggregate"]] <- Idents(object = P3F1)
  levels(P3F1) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')


# rename subtype
metadata <- P3F1@meta.data

# rename DS.Difference_common score into lineage score
colnames(P3F1@meta.data)[34] <- "Cluster_assignment"

# save metadata as df
metadata <- data.frame(P3F1@meta.data)

# load gene lists -----------------------------------
markers <- read.csv(file.path(genelist_dir, 'P3F1_cluster_markers.csv'))

# rename DNA replication as Proliferative
markers <- markers %>%
  mutate_all(~ ifelse(. == 'DNA replication', 'Proliferative', .))


# order by Annotation and fold change
markers <- markers %>% arrange(cluster, desc(avg_log2FC))
markers <- as.data.frame(markers)

# convert into list
marker_list <- split(markers[, -c(1:7)], f = markers$cluster)

# select top 10 genes per cluster
marker_list <- lapply(marker_list,head,50)


# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- marker_list
genes_df <- data.frame(
  Values = unlist(marker_list),
  Group = rep(names(marker_list), sapply(marker_list, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- P3F1@assays$RNA@data[genes_df$Values, ]


# define annotations -----------------------------------
# order dataframe based on clusters and then shuffle cells inside (random shuffling of Ds.difference)

order <- data.frame(metadata %>%
  arrange(Cluster_assignment, sample(row_number()))) 
ordered_rows <- rownames(order)

# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)

# order based on  score
seurat_df <- seurat_df[, ordered_rows]

# define annotations (important: order based on new order of cells)
metadata2 <- metadata[colnames(seurat_df), ]
metadata <- metadata2
annotation_info <- metadata[, c("subtype", "Cluster_assignment")]


# define colors  -----------------------------------
# Color Louvain clusters
col_Cluster_assignment <-  c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF' , '#D8D8E0FF')
names(col_Cluster_assignment) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')

col_subtype <- c('#D3A2C2FF', '#95CECFFF')
names(col_subtype) <- c('eRMS', 'aRMS')

col_column_ha <- list(subtype = col_subtype,
                      Cluster_assignment = col_Cluster_assignment)


# define column annotation
column_ha = HeatmapAnnotation(df = annotation_info,
                              col = col_column_ha)

# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)

# genes to mark
elements_to_find <- c('FN1', 'CAV1',  'MYOD1', 'MYOG', 'MYH3', 'MKI67', 'CENPF', 'TTN', 'MYL4', 'DCX', 'SYT1', 'L1CAM', 'STMN4')
rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1, 0, 1), c("#2E5A87FF", "#FCFDFEFF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_P3F1_heatmap_50.pdf"), width=10, height=6)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              #bottom_annotation = bottom_ha,
              border = TRUE,
              row_split = factor(c(genes_df$Group), 
                                 levels = c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')),
              top_annotation = column_ha,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun
)
draw(ht)
dev.off()






# -----------------------------------------------------------------------
# (3) Heatmap top marker genes FP-RMS (PAX7::FOXO1) clusters
# -----------------------------------------------------------------------

# load Seurat object -----------------------------------
P7F1 <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ARMS_P7F1_20230713.rds"))

# rename subtype
metadata <- P7F1@meta.data

# rename DS.Difference_common score into lineage score
colnames(P7F1@meta.data)[34] <- "Cluster_assignment"

# save metadata as df
metadata <- data.frame(P7F1@meta.data)

# load gene lists -----------------------------------
markers <- read.csv(file.path(genelist_dir, 'P7F1_cluster_markers.csv'))

# rename DNA replication as Proliferative
markers <- markers %>%
  mutate_all(~ ifelse(. == 'DNA replication', 'Proliferative', .))


# order by Annotation and fold change
markers <- markers %>% arrange(cluster, desc(avg_log2FC))
markers <- as.data.frame(markers)

# convert into list
marker_list <- split(markers[, -c(1:7)], f = markers$cluster)

# select top 10 genes per cluster
marker_list <- lapply(marker_list,head,50)


# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- marker_list
genes_df <- data.frame(
  Values = unlist(marker_list),
  Group = rep(names(marker_list), sapply(marker_list, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- P7F1@assays$RNA@data[genes_df$Values, ]


# define annotations -----------------------------------
# order dataframe based on clusters and then shuffle cells inside (random shuffling of Ds.difference)

order <- data.frame(metadata %>%
                      arrange(Cluster_assignment, sample(row_number()))) 
ordered_rows <- rownames(order)

# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)

# order based on  score
seurat_df <- seurat_df[, ordered_rows]

# define annotations (important: order based on new order of cells)
metadata2 <- metadata[colnames(seurat_df), ]
metadata <- metadata2
annotation_info <- metadata[, c("subtype", "Cluster_assignment")]


# define colors  -----------------------------------
# Color Louvain clusters
col_Cluster_assignment <-  c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#66C5CCFF')
names(col_Cluster_assignment) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal')

col_subtype <- c('#D3A2C2FF', '#95CECFFF')
names(col_subtype) <- c('eRMS', 'aRMS')

col_column_ha <- list(subtype = col_subtype,
                      Cluster_assignment = col_Cluster_assignment)


# define column annotation
column_ha = HeatmapAnnotation(df = annotation_info,
                              col = col_column_ha)

# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)

# genes to mark
elements_to_find <- c('FN1', 'CAV1',  'MYOD1', 'MYOG', 'MYH3', 'MKI67', 'CENPF', 'TTN', 'MYL4', 'DCX', 'SYT1', 'L1CAM', 'STMN4')
rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1, 0, 1), c("#2E5A87FF", "#FCFDFEFF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_P7F1_heatmap_50.pdf"), width=10, height=6)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              #bottom_annotation = bottom_ha,
              border = TRUE,
              row_split = factor(c(genes_df$Group), 
                                 levels = c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'Neuronal', 'Apoptosis')),
              top_annotation = column_ha,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun,
              row_names_rot = 30
)
draw(ht)
dev.off()










# -----------------------------------------------------------------------
# (4) Heatmap top marker genes FN-RMS clusters
# -----------------------------------------------------------------------

# load Seurat object -----------------------------------
ERMS <- readRDS(file.path(base_dir, "write/Danielli_Patel_Langenau_RPCA_ERMS_20230713.rds"))
  # rename ERMS (problem with current IDs)
  Idents(ERMS) = "cluster_names"
  new.cluster.ids.aggregate <- c('Progenitor', 'Progenitor', 'Proliferative',  'Differentiated', 'IFN',  'Proliferative', 'Ground', 'Ground')
  names(new.cluster.ids.aggregate) <- levels(ERMS)
  ERMS<- RenameIdents(ERMS, new.cluster.ids.aggregate)
  ERMS[["cluster_names_aggregate"]] <- Idents(object = ERMS)
  levels(ERMS) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'IFN')

# rename subtype
metadata <- ERMS@meta.data

# rename DS.Difference_common score into lineage score
colnames(ERMS@meta.data)[32] <- "Cluster_assignment"


# save metadata as df
metadata <- data.frame(ERMS@meta.data)

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

# select top 10 genes per cluster
marker_list <- lapply(marker_list,head,50)


# Subset seurat object to genes of interest -----------------------------------
# create df with genes of interest
genes <- marker_list
genes_df <- data.frame(
  Values = unlist(marker_list),
  Group = rep(names(marker_list), sapply(marker_list, length))
)
rownames(genes_df) <- NULL

# subset Seurat object to genes of interest 
seurat_df <- ERMS@assays$RNA@data[genes_df$Values, ]


# define annotations -----------------------------------
# order dataframe based on clusters and then shuffle cells inside (random shuffling of Ds.difference)

order <- data.frame(metadata %>%
                      arrange(Cluster_assignment, sample(row_number()))) 
ordered_rows <- rownames(order)

# scale dataset
seurat_df <- seurat_df - rowMeans(seurat_df)
seurat_df <- as.data.frame(seurat_df)

# order based on  score
seurat_df <- seurat_df[, ordered_rows]

# define annotations (important: order based on new order of cells)
metadata2 <- metadata[colnames(seurat_df), ]
metadata <- metadata2
annotation_info <- metadata[, c("subtype", "Cluster_assignment")]


# define colors  -----------------------------------
# Color Louvain clusters
col_Cluster_assignment <-  c("#7FBC41FF", '#FFAD72FF', '#FFE5CCFF', "#8E0152FF",'#B497E7FF')
names(col_Cluster_assignment) <- c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'IFN')

col_subtype <- c('#D3A2C2FF', '#95CECFFF')
names(col_subtype) <- c('eRMS', 'aRMS')

col_column_ha <- list(subtype = col_subtype,
                      Cluster_assignment = col_Cluster_assignment)


# define column annotation
column_ha = HeatmapAnnotation(df = annotation_info,
                              col = col_column_ha)

# define row annotation
row_ha = rowAnnotation(Group = genes_df$Group)

# genes to mark
elements_to_find <- c('FN1', 'CD44',  'MYOD1', 'MYOG', 'MYH3', 'MKI67', 'CENPF', 'TTN', 'MYL4', 'IFI6', 'SYT1', 'IFI44', 'ISG15')
rows_genes_to_mark <- which(genes_df$Values %in% elements_to_find)
rows_genes_to_mark_name <- genes_df$Values[rows_genes_to_mark]

right_ha = rowAnnotation(foo = anno_mark(at = rows_genes_to_mark, 
                                         labels = rows_genes_to_mark_name))

# define colors heatmap
col_fun = colorRamp2(c(-1, 0, 1), c("#2E5A87FF", "#FCFDFEFF", "#A90C38FF"))


# plot heatmap 
Cairo::CairoPDF(file.path(plot_dir, "1_ERMS_heatmap_50.pdf"), width=10, height=6)
ht <- Heatmap(seurat_df,
              cluster_columns = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE,
              #bottom_annotation = bottom_ha,
              border = TRUE,
              row_split = factor(c(genes_df$Group), 
                                 levels = c('Progenitor', 'Proliferative', 'Ground', 'Differentiated', 'IFN')),
              top_annotation = column_ha,
              right_annotation = right_ha,
              use_raster = TRUE,
              raster_quality = 10,
              cluster_rows = FALSE,
              col = col_fun,
              row_names_rot = 30
)
draw(ht)
dev.off()



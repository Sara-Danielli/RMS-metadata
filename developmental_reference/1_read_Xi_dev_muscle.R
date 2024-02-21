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
library(stringr)
options(Seurat.object.assay.version = "v5")

# Set up environment ------------------------------------------------
base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
Xi_dir <- file.path(base_dir, 'RMS/Xi_Cell_stem_cell_2020')
genelist_dir <- file.path(base_dir, 'list_final')

source(file.path('/Users/sdaniell/Dropbox (Partners HealthCare)/Sara Danielli/Manuscripts/2023 - Meta-data/GITHUB/RMS-metadata/Resources/Plot_style_v2.R'))

# color palette
#colors populations
col_celltype <- c('#C70E7BFF', '#FC6882FF',  '#6C6C9DFF',  '#A6E000FF','#1BB6AFFF','#172869FF')
names(col_celltype) <- c('Myogenic progenitors', 'Myoblasts-myocyte', 'Skeletal mesenchymal', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells')

#colors time points
col_timepoints <- c('#9E0142FF', '#D53E4FFF', '#F46D43FF', '#F46D43FF', '#FDAE61FF',  '#FEE08BFF',
                    '#FFFFBFFF',   '#E6F598FF', '#ABDDA4FF', '#66C2A5FF', '#3288BDFF', 
                    c(rep('#5E4FA2FF', 4)), 'grey')
names(col_timepoints) <- c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                           "Yr7",    "Yr11",   "Yr34" ,  "Yr42", NA)


# Load data ------------------------------------------------

# load human development datasets Xi et al. 2020 Stem Cell
week5_6 <- read.table(file.path(Xi_dir, 'myogenic_datasets/week5_6/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
week5_6_meta <- read.table(file.path(Xi_dir, 'myogenic_datasets/week5_6/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
week5_6 <- CreateSeuratObject(counts = week5_6, project = "Xi_2020")
week5_6 <- AddMetaData(week5_6, metadata = week5_6_meta)

week6_7 <- read.table(file.path(Xi_dir, 'myogenic_datasets/week6_7/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
week6_7_meta <- read.table(file.path(Xi_dir, 'myogenic_datasets/week6_7/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
week6_7 <- CreateSeuratObject(counts = week6_7, project = "Xi_2020")
week6_7 <- AddMetaData(week6_7, metadata = week6_7_meta)

week7_8 <- read.table(file.path(Xi_dir, 'myogenic_datasets/week7_8/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
week7_8_meta <- read.table(file.path(Xi_dir, 'myogenic_datasets/week7_8/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
week7_8 <- CreateSeuratObject(counts = week7_8, project = "Xi_2020")
week7_8 <- AddMetaData(week7_8, metadata = week7_8_meta)

week9 <- read.table(file.path(Xi_dir, 'myogenic_datasets/week9/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
week9_meta<- read.table(file.path(Xi_dir, 'myogenic_datasets/week9/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
week9 <- CreateSeuratObject(counts = week9, project = "Xi_2020")
week9 <- AddMetaData(week9, metadata = week9_meta)

week12_14 <- read.table(file.path(Xi_dir, 'myogenic_datasets/week12_14/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
week12_14_meta<- read.table(file.path(Xi_dir, 'myogenic_datasets/week12_14/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
week12_14 <- CreateSeuratObject(counts = week12_14, project = "Xi_2020")
week12_14 <- AddMetaData(week12_14, metadata = week12_14_meta)

week17_18 <- read.table(file.path(Xi_dir, 'myogenic_datasets/week17_18/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
week17_18_meta<- read.table(file.path(Xi_dir, 'myogenic_datasets/week17_18/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
week17_18 <- CreateSeuratObject(counts = week17_18, project = "Xi_2020")
week17_18 <- AddMetaData(week17_18, metadata = week17_18_meta)

juvenile <- read.table(file.path(Xi_dir, 'myogenic_datasets/juvenile/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
juvenile_meta <- read.table(file.path(Xi_dir, 'myogenic_datasets/juvenile/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
juvenile <- CreateSeuratObject(counts = juvenile, project = "Xi_2020")
juvenile  <- AddMetaData(juvenile , metadata = juvenile_meta)

adult <- read.table(file.path(Xi_dir, 'myogenic_datasets/adult/exprMatrix.tsv'), sep = '\t', header=TRUE, row.names=1)
adult_meta <- read.table(file.path(Xi_dir, 'myogenic_datasets/adult/meta.tsv'), sep = '\t', header=TRUE, row.names=1)
adult <- CreateSeuratObject(counts = adult, project = "Xi_2020")
adult <- AddMetaData(adult, metadata = adult_meta)


# merge data  ------------------------------------------------
merged_Xietal <- merge(week5_6, y = list(week6_7, week7_8, week9, week12_14, week17_18, juvenile, adult),
                       add.cell.ids = c('week5_6', 'week6_7', 'week7_8', 'week9', 'week12_14', 'week17_18', 'juvenile', 'adult'),
                       project = "development_Xi_2020")
merged_Xietal <- JoinLayers(merged_Xietal)

Idents(merged_Xietal) <- 'orig.ident.short'

DefaultAssay(merged_Xietal) <- "RNA"

# Define mitochondrial genes
merged_Xietal[["percent.mt"]] <- PercentageFeatureSet(merged_Xietal, pattern = "^MT-")

# Refine metadata annotations  ------------------------------------------------
# (1) add cell type info
cellname <- colnames(merged_Xietal)
celltype <- unlist(str_split(merged_Xietal@meta.data$Cell.Type, "-", n = 1))
#create dataframe
df <- data.frame(cellname, celltype)
# add info as metadata  
merged_Xietal <- AddMetaData(merged_Xietal, metadata = df)

# add average cell type (MP, MB, ecc) and time point
celltype <- c()
timepoint <- c()
# extract average cell type
cell_type <- merged_Xietal@meta.data$orig.ident.short_cell.type
# Extract the third element
for (i in 1:length(cell_type)) {
  split_elements <- str_split(cell_type[i], "_", n = 3)
  celltype[i] <- split_elements[[1]][3]
  timepoint[i] <- split_elements[[1]][1]
}
#create dataframe
df <- data.frame(cellname, celltype, timepoint)
# add info as metadata  
merged_Xietal <- AddMetaData(merged_Xietal, metadata = df)


# Rename celltype names  -----------------------------------
metadata <- merged_Xietal@meta.data
# Use dplyr's mutate to substitute specific cells with new names
metadata <- metadata %>%
  mutate(celltype = case_when(
    celltype == "MP" ~ "Myogenic progenitors",
    celltype == "MB-MC" ~ "Myoblasts-myocyte",
    celltype == "SkM.Mesen" ~ "Skeletal mesenchymal",
    celltype == "MB" ~ "Myoblasts",
    celltype == "MC" ~ "Myocytes",
    celltype == "SC" ~ "Postnatal satellite cells"
  ))

# Add new metadata to Seurat objetct
merged_Xietal <- AddMetaData(merged_Xietal, metadata)

# Downstream analysis Seurat object ------------------------------------------------

## Data normalization
merged_Xietal <- NormalizeData(merged_Xietal, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable features
merged_Xietal <- FindVariableFeatures(merged_Xietal, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(merged_Xietal), 10)

## Load signatures
Common_differentiated <- read_excel(file.path(genelist_dir, "Common_differentiated.xlsx"), col_names=FALSE)
Common_differentiated <- as.list(Common_differentiated)
Common_proliferative <- read_excel(file.path(genelist_dir,"Common_proliferative.xlsx"), col_names=FALSE)
Common_proliferative <- as.list(Common_proliferative)
Common_Stemcell <- read_excel(file.path(genelist_dir, "Common_Stemcell.xlsx"), col_names=FALSE)
Common_Stemcell <- as.list(Common_Stemcell)

signatures <- c(Common_differentiated, Common_proliferative, Common_Stemcell)
names(signatures) <- list('Differentiated score', 'Proliferative score', 'Progenitor score')

# select top 25 genes per signature
signatures <- lapply(signatures,head,25)

# Score for metaprograms  ----------------------------------
merged_Xietal <- AddModuleScore(object = merged_Xietal, assay = 'RNA', features = signatures, name = names(signatures))
merged_Xietal <- ScaleData(merged_Xietal)

# rename metadata names of scores
col_start <- length(colnames(merged_Xietal@meta.data)) - length(names(signatures)) + 1
# identify number of last column with metadata scores
col_end <- length(colnames(merged_Xietal@meta.data))
# rename columns with score name
colnames(merged_Xietal@meta.data)[col_start:col_end] <- names(signatures)


## Dimension reduction
merged_Xietal
merged_Xietal <- RunPCA(merged_Xietal, reduction.name = "pca_no_regression", reduction.key = "pca_no_regression")

# reorder labels
merged_Xietal$timepoint <- factor(x = merged_Xietal$timepoint, 
                                  levels = c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                                             "Yr7",    "Yr11",   "Yr34" ,  "Yr42"))

merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                 levels = c("Myogenic progenitors", "Myoblasts-myocyte", "Myoblasts",  "Myocytes",
                                            "Skeletal mesenchymal",  "Postnatal satellite cells")) 

# Plot PCA
DimPlot(merged_Xietal, reduction = "pca_no_regression", group.by = 'celltype', cols = col_celltype) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/1_PCA_celltype.pdf", 
       width=10, height=5, dpi=300)

ElbowPlot(merged_Xietal, reduction='pca_no_regression')
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/2_Elbowplot.pdf", 
       width=6, height=4, dpi=300)


# Cluster cells
merged_Xietal <- FindNeighbors(merged_Xietal, reduction = 'pca_no_regression', dims = 1:10)

# Plot UMAP
merged_Xietal <- RunUMAP(merged_Xietal, dims = 1:10, reduction = 'pca_no_regression', reduction.name = 'UMAP_no_regression')

# save object  ------------------------------------------------------------
saveRDS(merged_Xietal, file = "/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")

merged_Xietal <- readRDS("/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")



## Plots ------------------------------------------------------------
DimPlot(merged_Xietal, reduction = "UMAP_no_regression", group.by = 'timepoint', cols = col_timepoints) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/3_UMAP_celltype_col1.pdf", 
       width=12, height=5, dpi=300)

DimPlot(merged_Xietal, reduction = "UMAP_no_regression", group.by = 'celltype', cols=col_celltype) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/3_UMAP_celltype_col2.pdf", 
       width=12, height=5, dpi=300)


# UMAP plot of RMS scores
p <- FeaturePlot(merged_Xietal, features = names(signatures), reduction = "UMAP_no_regression", combine = FALSE, order = TRUE, pt.size = 2, raster=TRUE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoAxes() + scale_colour_distiller(palette="RdBu")
}
cowplot::plot_grid(plotlist = p, ncol=3)
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/4_UMAP_scores.pdf", 
       width=10, height=3, dpi=300)

# Dotplot of RMS scores
merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                 levels = c("Postnatal satellite cells",  "Skeletal mesenchymal", "Myocytes", "Myoblasts",
                                            "Myoblasts-myocyte", "Myogenic progenitors"))

DotPlot(merged_Xietal, features = rev(names(signatures)),
        group.by = 'celltype',
        assay = 'RNA', 
        scale.min = 100,
        scale.max = 100,
        dot.scale=10,
        cols = c("white", "red3"),
        scale = FALSE)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/5_DotPlot_scores.pdf", 
       width=5, height=3, dpi=300)

DotPlot(merged_Xietal, 
        group.by = 'timepoint',
        features = rev(names(signatures)),
        assay = 'RNA', 
        scale.min = 100,
        scale.max = 100,
        dot.scale=10,
        cols = c("white", "red3"),
        scale = FALSE)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/5_DotPlot_scores_timepoints.pdf", 
       width=4, height=5, dpi=300)



# Identify cluster markers (by cell type and time point) ------------------------------------------

# define celltype markers and select top 25 marker genes
merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                 levels = c("Myogenic progenitors", "Myoblasts-myocyte", "Myoblasts",  "Myocytes",
                                            "Skeletal mesenchymal",  "Postnatal satellite cells")) 
Idents(merged_Xietal) = 'celltype'
muscle.markers <- FindAllMarkers(merged_Xietal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3, test.use = 'negbinom')

genes_per_cluster = 25
top_genes <- muscle.markers %>%
  group_by(cluster) %>%
  slice_max(n = genes_per_cluster, order_by = avg_log2FC)
# save as csv
write.csv(muscle.markers, "/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/7_markers_celltype.csv")
# save as list
muscle.markers <- muscle.markers %>% arrange(cluster, desc(avg_log2FC))
muscle.markers_list <- split(muscle.markers[, -c(1:6)], f = muscle.markers$cluster)
saveRDS(muscle.markers_list, file.path(base_dir, 'output/metadata/dev_muscle/7_markers_celltype_list.rds'))


# define celltime point  markers and select top 25 marker genes
Idents(merged_Xietal) <- 'timepoint'
muscle.time.markers <- FindAllMarkers(merged_Xietal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3, test.use = 'negbinom')

genes_per_cluster = 25
top_genes.cell <- muscle.time.markers %>%
  group_by(cluster) %>%
  slice_max(n = genes_per_cluster, order_by = avg_log2FC)
# save as csv
write.csv(muscle.time.markers, "/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/7_markers_timepoint.csv")
# save as list
muscle.time.markers <- muscle.time.markers %>% arrange(cluster, desc(avg_log2FC))
muscle.time.markers_list <- split(muscle.time.markers[, -c(1:6)], f = muscle.time.markers$cluster)
saveRDS(muscle.time.markers_list, file.path(base_dir, 'output/metadata/dev_muscle/7_markers_timepoint_list.rds'))


# Plot cell type distribution ------------------------------------------

# reorder labels
merged_Xietal$timepoint <- factor(x = merged_Xietal$timepoint, 
                                  levels = c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                                             "Yr7",    "Yr11",   "Yr34" ,  "Yr42"))

merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                 levels = c("Myogenic progenitors", "Myoblasts-myocyte", "Myoblasts",  "Myocytes",
                                            "Skeletal mesenchymal",  "Postnatal satellite cells")) 


plot_bar(merged_Xietal, merged_Xietal$timepoint, merged_Xietal$celltype, col = col_celltype) 
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/4_cell_composition_pop.pdf", 
       width=8, height=4)

plot_bar(merged_Xietal, merged_Xietal$celltype, merged_Xietal$timepoint, col = col_timepoints) 
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/4_cell_composition_timepoints.pdf", 
       width=4, height=5)


# Scoring cells for FN-RMS metaprograms   ------------------------------------------
merged_Xietal <- readRDS("/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")

## Data normalization
merged_Xietal <- NormalizeData(merged_Xietal, normalization.method = "LogNormalize", scale.factor = 10000)
## Highly variable features
merged_Xietal <- FindVariableFeatures(merged_Xietal, selection.method = "vst", nfeatures = 2000)

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

# select top 25 genes per cluster
marker_list <- lapply(marker_list,head,25)


# Score for metaprograms  --------------------
merged_Xietal <- AddModuleScore(object = merged_Xietal, assay = 'RNA', features = marker_list, name = names(marker_list))
merged_Xietal <- ScaleData(merged_Xietal)

# rename metadata names of scores
col_start <- length(colnames(merged_Xietal@meta.data)) - length(names(marker_list)) + 1
# identify number of last column with metadata scores
col_end <- length(colnames(merged_Xietal@meta.data))
# rename columns with score name
colnames(merged_Xietal@meta.data)[col_start:col_end] <- names(marker_list)

merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                     levels = c("Myogenic progenitors", "Myoblasts-myocyte", "Myoblasts",                
                                                "Myocytes", "Skeletal mesenchymal", "Postnatal satellite cells"))

DotPlot(merged_Xietal, features = c('IFN', 'Differentiated', 'Ground', 'Proliferative', 'Progenitor'),
        group.by = 'celltype',
        assay = 'RNA', 
        scale.min = 100,
        scale.max = 100,
        dot.scale=10,
        cols = c("white", "red3"),
        scale = FALSE)  + 
  #scale_colour_distiller(palette="RdBu") +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1, colour="black"),
        axis.text.y = element_text(size=12, colour="black"),
        axis.title=element_blank(),
        legend.text = element_text(size = 12),
        #legend.title = element_blank()
  ) 
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/6_DotPlot_scores_FNRMS.pdf", 
       width=5, height=3, dpi=300)



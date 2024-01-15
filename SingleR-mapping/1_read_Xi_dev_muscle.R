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

##################################################################
#  (0) Set up directories and environment:  #
##################################################################

base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
Xi_dir <- file.path(base_dir, 'RMS/Xi_Cell_stem_cell_2020')
genelist_dir <- file.path(base_dir, 'list_final')

# color scheme
#colors populations
col_celltype <- c('#C70E7BFF', '#FC6882FF',  '#6C6C9DFF',  '#A6E000FF','#1BB6AFFF','#172869FF')
names(col_celltype) <- c('Myogenic progenitors', 'Myoblasts-myocyte', 'Skeletal mesenchymal', 'Myoblasts', 'Myocytes', 'Postnatal satellite cells')

#colors time points
col_timepoints <- c('#9E0142FF', '#D53E4FFF', '#F46D43FF', '#F46D43FF', '#FDAE61FF',  '#FEE08BFF',
                     '#FFFFBFFF',   '#E6F598FF', '#ABDDA4FF', '#66C2A5FF', '#3288BDFF', 
                    c(rep('#5E4FA2FF', 4)))
names(col_timepoints) <- c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                           "Yr7",    "Yr11",   "Yr34," ,  "Yr42")


##################################################################
#  (1) Load human muscle dataset Xi et al. 2020:  #
##################################################################

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



merged_Xietal <- merge(week5_6, y = list(week6_7, week7_8, week9, week12_14, week17_18, juvenile, adult),
                       add.cell.ids = c('week5_6', 'week6_7', 'week7_8', 'week9', 'week12_14', 'week17_18', 'juvenile', 'adult'),
                       project = "development_Xi_2020")
merged_Xietal <- JoinLayers(merged_Xietal)

Idents(merged_Xietal) <- 'orig.ident.short'

DefaultAssay(merged_Xietal) <- "RNA"

# Define mitochondrial genes
merged_Xietal[["percent.mt"]] <- PercentageFeatureSet(merged_Xietal, pattern = "^MT-")

# Refine metadata annotations
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
  
#save object  
saveRDS(merged_Xietal, file = "/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")
  
merged_Xietal <- readRDS("/Volumes/Sara_PhD/scRNAseq_data/write/Dev_muscle_Xietal.rds")

##################################################################
#  (2) Downstream analyses:  #
##################################################################

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

signatures <- list(Common_differentiated, Common_proliferative, Common_Stemcell)
names(signatures) <- list('RMS_differentiated', 'RMS_proliferative', 'RMS_progenitor')

## score
merged_Xietal <- AddModuleScore(object = merged_Xietal, features = Common_differentiated, name = 'Differentiated')
merged_Xietal <- AddModuleScore(object = merged_Xietal, features = Common_proliferative, name = 'Proliferative')
merged_Xietal <- AddModuleScore(object = merged_Xietal, features = Common_Stemcell, name = 'Progenitor')


# Scaling data
merged_Xietal <- ScaleData(merged_Xietal)

## Dimension reduction
merged_Xietal
merged_Xietal <- RunPCA(merged_Xietal, reduction.name = "pca_no_regression", reduction.key = "pca_no_regression")

# reorder labels
merged_Xietal$timepoint <- factor(x = merged_Xietal$timepoint, 
                               levels = c("Wk5.0", "Wk6.0" ,  "Wk6.5", "Wk7.0", "Wk7.25",  "Wk7.75" , "Wk9" ,   "Wk12",   "Wk14", "Wk17" ,  "Wk18", 
                                          "Yr7",    "Yr11",   "Yr34," ,  "Yr42"))

merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                 levels = c("Myogenic progenitors", "Myoblasts-myocyte", "Myoblasts",  "Myocytes",
                                            "Skeletal mesenchymal",  "Postnatal satellite cells")) 

# PCA plots
DimPlot(merged_Xietal, reduction = "pca_no_regression", group.by = 'celltype', cols = col_celltype) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/1_PCA_celltype.pdf", 
       width=10, height=5, dpi=300)

ElbowPlot(merged_Xietal, reduction='pca_no_regression')
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/2_Elbowplot.pdf", 
       width=6, height=4, dpi=300)


## Cluster cells
merged_Xietal <- FindNeighbors(merged_Xietal, reduction = 'pca_no_regression', dims = 1:10)

# UMAP
merged_Xietal <- RunUMAP(merged_Xietal, dims = 1:10, reduction = 'pca_no_regression', reduction.name = 'UMAP_no_regression')

## Plots
DimPlot(merged_Xietal, reduction = "UMAP_no_regression", group.by = 'timepoint', cols = col_timepoints) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/3_UMAP_celltype_col1.pdf", 
       width=12, height=5, dpi=300)

DimPlot(merged_Xietal, reduction = "UMAP_no_regression", group.by = 'celltype', cols=col_celltype) + NoAxes()
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/3_UMAP_celltype_col2.pdf", 
       width=12, height=5, dpi=300)


# UMAP plot of RMS scores
scores <- c('Progenitor1', 'Proliferative1', 'Differentiated1')
p <- FeaturePlot(merged_Xietal, features = scores, reduction = "UMAP_no_regression", combine = FALSE, order = TRUE, pt.size = 2, raster=TRUE)
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
Idents(merged_Xietal) = 'celltype'

DotPlot(merged_Xietal, features = scores, 
        assay = 'RNA', 
        col.min = 0, 
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
       width=5, height=4, dpi=300)



########### FIND CLUSTER MARKERS ###########
merged_Xietal$celltype <- factor(x = merged_Xietal$celltype, 
                                 levels = c("Myogenic progenitors", "Myoblasts-myocyte", "Myoblasts",  "Myocytes",
                                            "Skeletal mesenchymal",  "Postnatal satellite cells")) 
Idents(merged_Xietal) = 'celltype'


merged_Xietal.markers <- FindAllMarkers(merged_Xietal, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 1)
write.csv(merged_Xietal.markers, "/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/7_markers.csv")

# save as list
marker_list <- split(merged_Xietal.markers[, -c(1:6)], f = merged_Xietal.markers$cluster)
saveRDS(marker_list, file.path(base_dir, 'output/metadata/dev_muscle/7_markers.rds'))


plot_bar <-ggplot(merged_Xietal@meta.data, aes(x=orig.ident.short, fill= (merged_Xietal@meta.data$celltype))) + 
  geom_bar(position = "fill", color="black") +
  scale_fill_manual(values = col_celltype) +
  labs (y='Proportion', x='') + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot_bar
ggsave("/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/4_cell_composition_pop.pdf", 
       width=10, height=4, dpi=300)

## Calculate number of cells within each cluster
library(data.table)
md <- merged_Xietal@meta.data %>% as.data.table
md[, .N, by = c("orig.ident.short", "orig.ident.short_cell.type")]
table_PDX <- md[, .N, by = c("orig.ident.short", "orig.ident.short_cell.type")] %>% dcast(., orig.ident.short ~ orig.ident.short_cell.type, value.var = "N")
write.csv(table_PDX, "/Volumes/Sara_PhD/scRNAseq_data/output/metadata/dev_muscle/5_Cluster_information_name.csv")







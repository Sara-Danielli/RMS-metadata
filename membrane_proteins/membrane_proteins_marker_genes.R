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
library(writexl)
library(stringi)
library(stringr)


# Organize environment  -----------------------------------
base_dir <- '/Volumes/Sara_PhD/scRNAseq_data'

genelist_dir <- file.path(base_dir, 'list_final')

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/membrane_proteins')

plot_dir <- file.path(analysis_dir, 'plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
write_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}


# load list of known plasma proteins -----------------------------------
proteins <- read.csv(file.path(genelist_dir, 'plasma_membrane_proteins/subcellular_location.csv'))

# filter for plasma membrane (try to get all that contain membrane)
#membrane_proteins <- proteins %>% filter(Main.location =='Plasma membrane')
membrane_proteins <- proteins %>%
  filter(str_detect(Main.location, 'membrane'))


membrane_proteins_genes <- membrane_proteins$Gene.name

# -----------------------------------------------------------------------
# (1) RMS atlas
# -----------------------------------------------------------------------

# load gene lists -----------------------------------
markers <- read_xlsx(file.path(genelist_dir, 'all_integrated_cluster_markers.xlsx'))
# order by Annotation and fold change
markers <- markers %>% arrange(Annotation, 'Average log2 fold change')
markers <- as.data.frame(markers)

# convert into list
marker_list <- split(markers[, -c(2:7)], f = markers$Annotation)

# Remove duplicates within each vector
marker_list <- lapply(marker_list, unique)

# select top 100 genes
#marker_list <- lapply(marker_list,head,100)

count_total <- c()
count_membrane_genes <- c()
name_membrane_genes  <- c()

for (i in seq_along(names(marker_list))) {
  count_total[i] <- length(marker_list[[i]])
  count_membrane_genes[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  name_membrane_genes[[i]] <- as.vector(marker_list[[i]][marker_list[[i]] %in% membrane_proteins_genes])
  #frequency_membrane_genes[i] <- count_membrane_genes[i]/length(marker_list[[i]])*100
}
names(count_total) <- names(count_membrane_genes) <- names(name_membrane_genes)  <-names(marker_list) #<- names(frequency_membrane_genes) 
name_membrane_genes <- stri_list2matrix(name_membrane_genes)
colnames(name_membrane_genes) <- names(marker_list)


# -----------------------------------------------------------------------
# (2) Heatmap top marker genes FP-RMS (PAX3::FOXO1) clusters
# -----------------------------------------------------------------------

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

# Remove duplicates within each vector
marker_list <- lapply(marker_list, unique)

# select top 100 genes
#marker_list <- lapply(marker_list,head,100)

count_total_P3F1 <- c()
count_membrane_genes_P3F1 <- c()
name_membrane_genes_P3F1 <- c()

for (i in seq_along(names(marker_list))) {
  count_total_P3F1[i] <- length(marker_list[[i]])
  count_membrane_genes_P3F1[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  name_membrane_genes_P3F1[[i]] <- as.vector(marker_list[[i]][marker_list[[i]] %in% membrane_proteins_genes])
  
}
names(count_total_P3F1) <- names(count_membrane_genes_P3F1) <- names(name_membrane_genes_P3F1) <- names(marker_list)
name_membrane_genes_P3F1 <- stri_list2matrix(name_membrane_genes_P3F1)
colnames(name_membrane_genes_P3F1) <- names(marker_list)

# -----------------------------------------------------------------------
# (3) Heatmap top marker genes FP-RMS (PAX7::FOXO1) clusters
# -----------------------------------------------------------------------

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

# Remove duplicates within each vector
marker_list <- lapply(marker_list, unique)

# select top 100 genes
#marker_list <- lapply(marker_list,head,100)

count_total_P7F1 <- c()
count_membrane_genes_P7F1 <- c()
name_membrane_genes_P7F1 <- c()

for (i in seq_along(names(marker_list))) {
  count_total_P7F1[i] <- length(marker_list[[i]])
  count_membrane_genes_P7F1[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  name_membrane_genes_P7F1[[i]] <- as.vector(marker_list[[i]][marker_list[[i]] %in% membrane_proteins_genes])
  
}
names(count_total_P7F1) <- names(count_membrane_genes_P7F1) <- names(name_membrane_genes_P7F1) <- names(marker_list)

name_membrane_genes_P7F1 <- stri_list2matrix(name_membrane_genes_P7F1)
colnames(name_membrane_genes_P7F1) <- names(marker_list)

# -----------------------------------------------------------------------
# (4) Heatmap top marker genes FN-RMS clusters
# -----------------------------------------------------------------------

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

# select top 100 genes
#marker_list <- lapply(marker_list,head,100)

count_total_ERMS <- c()
count_membrane_genes_ERMS <- c()
name_membrane_genes_ERMS <- c()

for (i in seq_along(names(marker_list))) {
  count_total_ERMS[i] <- length(marker_list[[i]])
  count_membrane_genes_ERMS[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  name_membrane_genes_ERMS[[i]] <- as.vector(marker_list[[i]][marker_list[[i]] %in% membrane_proteins_genes])
}
names(count_total_ERMS) <- names(count_membrane_genes_ERMS) <- names(name_membrane_genes_ERMS) <- names(marker_list)

name_membrane_genes_ERMS <- stri_list2matrix(name_membrane_genes_ERMS)
colnames(name_membrane_genes_ERMS) <- names(marker_list)

# -----------------------------------------------------------------------
# Put together results
# -----------------------------------------------------------------------

counts_membrane <- list(count_membrane_genes, count_membrane_genes_P3F1, count_membrane_genes_P7F1, count_membrane_genes_ERMS)
name_membrane_genes_combined <- list(name_membrane_genes, name_membrane_genes_P3F1, name_membrane_genes_P7F1, name_membrane_genes_ERMS)
names(counts_membrane) <- names(name_membrane_genes_combined) <- c('RMS atlas', 'PAX3-FOXO1', 'PAX7-FOXO1', 'FN-RMS')

# export name of genes encoding for membrane proteins
df <- list()
name_membrane_genes_combined2 <- list()
for (i in seq_along(counts_membrane)){
df[[i]] <- cbind(as.data.frame(counts_membrane[[i]]))
name_membrane_genes_combined2[[i]] <- as.data.frame(name_membrane_genes_combined[[i]])
}
names(name_membrane_genes_combined2) <- names(df) <-  c('RMS atlas', 'PAX3-FOXO1', 'PAX7-FOXO1', 'FN-RMS')
write_xlsx(name_membrane_genes_combined2, file.path(write_dir, paste0('names_membrane_genes.xlsx')))

# export numbers of genes encoding for membraneproteins
for (i in seq_along(df)){
  write.csv(df[[i]], file.path(write_dir, paste0('number_membrane_genes', names(df)[i], '.csv')))
}

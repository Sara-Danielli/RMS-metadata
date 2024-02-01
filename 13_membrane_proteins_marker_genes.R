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

analysis_dir <- file.path(base_dir, 'output/metadata/Patel_Danielli_Langenau/RPCA_name/membrane_proteins')

plot_dir <- file.path(analysis_dir, 'plot')
if (!dir.exists(plot_dir)){dir.create(plot_dir, recursive = T)}
write_dir <- file.path(analysis_dir, 'data')
if (!dir.exists(write_dir)){dir.create(write_dir, recursive = T)}


# load list of known plasma proteins -----------------------------------
proteins <- read_csv(file.path(genelist_dir, 'plasma_membrane_proteins/subcellular_location.csv'))

# filter for plasma membrane
membrane_proteins <- proteins %>% filter(`Main location` =='Plasma membrane')
membrane_proteins_genes <- membrane_proteins$`Gene name`

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

count_total <- c()
count_membrane_genes <- c()
frequency_membrane_genes <- c()

for (i in seq_along(names(marker_list))) {
  count_total[i] <- length(marker_list[[i]])
  count_membrane_genes[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  frequency_membrane_genes[i] <- count_membrane_genes[i]/length(marker_list[[i]])*100
}
names(count_total) <- names(count_membrane_genes) <- names(frequency_membrane_genes) <-names(marker_list)


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

count_total_P3F1 <- c()
count_membrane_genes_P3F1 <- c()
frequency_membrane_genes_P3F1 <- c()

for (i in seq_along(names(marker_list))) {
  count_total_P3F1[i] <- length(marker_list[[i]])
  count_membrane_genes_P3F1[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  frequency_membrane_genes_P3F1[i] <- count_membrane_genes_P3F1[i]/length(marker_list[[i]])*100
}
names(count_total_P3F1) <- names(count_membrane_genes_P3F1) <- names(frequency_membrane_genes_P3F1) <- names(marker_list)


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

count_total_P7F1 <- c()
count_membrane_genes_P7F1 <- c()
frequency_membrane_genes_P7F1 <- c()

for (i in seq_along(names(marker_list))) {
  count_total_P7F1[i] <- length(marker_list[[i]])
  count_membrane_genes_P7F1[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  frequency_membrane_genes_P7F1[i] <- count_membrane_genes_P7F1[i]/length(marker_list[[i]])*100
}
names(count_total_P7F1) <- names(count_membrane_genes_P7F1) <- names(frequency_membrane_genes_P7F1) <- names(marker_list)



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

count_total_ERMS <- c()
count_membrane_genes_ERMS <- c()
frequency_membrane_genes_ERMS <- c()

for (i in seq_along(names(marker_list))) {
  count_total_ERMS[i] <- length(marker_list[[i]])
  count_membrane_genes_ERMS[i] <- sum(marker_list[[i]] %in% membrane_proteins_genes)
  frequency_membrane_genes_ERMS[i] <- count_membrane_genes_ERMS[i]/length(marker_list[[i]])*100
}
names(count_total_ERMS) <- names(count_membrane_genes_ERMS) <- names(frequency_membrane_genes_ERMS) <- names(marker_list)




# -----------------------------------------------------------------------
# Put together results
# -----------------------------------------------------------------------

counts_membrane <- list(count_membrane_genes, count_membrane_genes_P3F1, count_membrane_genes_P7F1, count_membrane_genes_ERMS)
counts_total <- list(count_total, count_total_P3F1, count_total_P7F1, count_total_ERMS)
frequency <- list(frequency_membrane_genes, frequency_membrane_genes_P3F1, frequency_membrane_genes_P7F1, frequency_membrane_genes_ERMS)

names(frequency) <- names(counts_membrane) <- names(counts_total) <- c('RMS atlas', 'PAX3::FOXO1', 'PAX7::FOXO1', 'FN-RMS')

# bind into df
df <- list()
for (i in seq_along(counts_membrane)){
df[[i]] <- cbind(as.data.frame(counts_membrane[[i]]), as.data.frame(counts_total[[i]]), as.data.frame(frequency[[i]]))
}

names(df) <- names(frequency) <- names(counts_membrane) <- names(counts_total) <- c('RMS atlas', 'PAX3::FOXO1', 'PAX7::FOXO1', 'FN-RMS')

for (i in seq_along(counts_membrane)){
  df[[i]] <- t(df[[i]])
  barplot(df, main = title, names.arg = seq_along(data_vector), col = "skyblue", ylim = c(0, max(data_vector) + 5))
  text(seq_along(data_vector), data_vector + 1, labels_vector, pos = 3, col = "red")
}


rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(data.table)
library(SeuratObject)
library(pheatmap)
library(ComplexHeatmap)

base_dir <- file.path('/Volumes/Sara_PhD/scRNAseq_data')
genelist_dir <- file.path(base_dir, 'list_final')

# Upsetplots
genes <- read_xlsx(file.path(genelist_dir, "UpSetPlot_progenitor_input.xlsx"))
genes <- genes[3:5]
genes_list <- lapply(genes, na.omit)

# convert to matrix
genes_matrix <- list_to_matrix(genes_list)

# make comb mat
#m_progenitor = make_comb_mat(genes_matrix, mode = 'intersect')
m_progenitor = make_comb_mat(genes_matrix, mode = 'distinct')

# make plot
pdf(file.path(genelist_dir, "progenitor.pdf"), width=4, height=4)
ht <- UpSet(m_progenitor[comb_degree(m_progenitor) >= 2])
draw(ht)
dev.off()

# Upsetplots 2
genes <- read_xlsx(file.path(genelist_dir, "UpSetPlot_proliferative_input.xlsx"))
genes_list <- lapply(genes, na.omit)

# convert to matrix
genes_matrix <- list_to_matrix(genes_list)

# make comb mat
#m_proliferative = make_comb_mat(genes_matrix, mode = 'intersect')
m_proliferative = make_comb_mat(genes_matrix, mode = 'distinct')

# make plot
pdf(file.path(genelist_dir, "proliferative.pdf"), width=4, height=4)
ht <- UpSet(m_proliferative[comb_degree(m_proliferative) >= 2])
draw(ht)
dev.off()

# Upsetplots 3
genes <- read_xlsx(file.path(genelist_dir, "UpSetPlot_differentiated_input.xlsx"))
genes_list <- lapply(genes, na.omit)

# convert to matrix
genes_matrix <- list_to_matrix(genes_list)

# make comb mat
#m_differentiated = make_comb_mat(genes_matrix, mode = 'intersect')
m_differentiated = make_comb_mat(genes_matrix, mode = 'distinct')

# make plot
pdf(file.path(genelist_dir, "differentiated.pdf"), width=4, height=4)
ht <- UpSet(m_differentiated[comb_degree(m_differentiated) >= 2])
draw(ht)
dev.off()


# combine 
top_ha = HeatmapAnnotation(
  "Progenitor" = anno_barplot(comb_size(m_progenitor), 
                           gp = gpar(fill = "black"), add_numbers = TRUE, height = unit(3, "cm")), 
  "Proliferative" = anno_barplot(comb_size(m_proliferative), 
                             gp = gpar(fill = "black"), add_numbers = TRUE, height = unit(3, "cm")), 
  "Differentiated" = anno_barplot(comb_size(m_differentiated), 
                         gp = gpar(fill = "black"), add_numbers = TRUE, height = unit(3, "cm")), 
  gap = unit(2, "mm"), annotation_name_side = "left", annotation_name_rot = 0)


# make plot
pdf(file.path(genelist_dir, "combined.pdf"), width=6, height=6)
ht <- UpSet(m_progenitor, top_annotation = top_ha, pt_size = unit(5, "mm"))
draw(ht)
dev.off()





## Subset subtypes
scRNAseq <- subset(PDX.integrated, subset = sequencing == 'scRNAseq')

scRNAseq <- ScaleData(scRNAseq, verbose = TRUE)

scores <- list(scRNAseq$Wei_mesenchymal1, scRNAseq$Wei_proliferative1, scRNAseq$Wei_muscle1,
               scRNAseq$Danielli_MuSC1,  scRNAseq$Danielli_cycling1, scRNAseq$Danielli_differentiated1,
               scRNAseq$Patel_mesoderm1, scRNAseq$Patel_myoblast1, scRNAseq$Patel_myocyte1,
               scRNAseq$Common_Stemcell1, scRNAseq$Common_proliferative1, scRNAseq$Common_differentiated1)


cor_res <- list()

for (i in 1:length(scores)) {
  cor_res[i] <- cor(scores[[i]], scRNAseq$Common_Stemcell1)
}

names(cor_res) <- c('Wei_mesenchymal1', "Wei_proliferative1", "Wei_muscle1",
                    "Danielli_MuSC1",  "Danielli_cycling1", "Danielli_differentiated1",
                    "Patel_mesoderm1", "Patel_myoblast1", "Patel_myocyte1",
                    "Common_Stemcell1", "Common_proliferative1", "Common_differentiated1")
cor_res <- as.data.frame(cor_res)
print(cor_res)
write.csv(file.path(analysis_dir, 'Correlation_common_progenitor.csv'))

cor_res <- list()

for (i in 1:length(scores)) {
  cor_res[i] <- cor(scores[[i]], scRNAseq$Common_proliferative1)
}

names(cor_res) <- c('Wei_mesenchymal1', "Wei_proliferative1", "Wei_muscle1",
                    "Danielli_MuSC1",  "Danielli_cycling1", "Danielli_differentiated1",
                    "Patel_mesoderm1", "Patel_myoblast1", "Patel_myocyte1",
                    "Common_Stemcell1", "Common_proliferative1", "Common_differentiated1")
cor_res <- as.data.frame(cor_res)
print(cor_res)
write.csv(file.path(analysis_dir, 'Correlation_common_proliferative.csv'))


cor_res <- list()

for (i in 1:length(scores)) {
  cor_res[i] <- cor(scores[[i]], scRNAseq$Common_proliferative1)
}

names(cor_res) <- c('Wei_mesenchymal1', "Wei_proliferative1", "Wei_muscle1",
                    "Danielli_MuSC1",  "Danielli_cycling1", "Danielli_differentiated1",
                    "Patel_mesoderm1", "Patel_myoblast1", "Patel_myocyte1",
                    "Common_Stemcell1", "Common_proliferative1", "Common_differentiated1")
cor_res <- as.data.frame(cor_res)
print(cor_res)
write.csv(file.path(analysis_dir, 'Correlation_common_proliferative.csv'))








## Publication-specific scores
Wei_mesenchymal <- read_excel(file.path(genelist_dir, "Wei_2022/Mesenchymal.xlsx"), col_names=FALSE)
Wei_mesenchymal <- unlist(Wei_mesenchymal)
Wei_proliferative <- read_excel(file.path(genelist_dir, "Wei_2022/Proliferative.xlsx"), col_names=FALSE)
Wei_proliferative <- unlist(Wei_proliferative)
Wei_muscle <- read_excel(file.path(genelist_dir, "Wei_2022/Muscle.xlsx"), col_names=FALSE)
Wei_muscle <- unlist(Wei_muscle)

Danielli_MuSC <- read_excel(file.path(genelist_dir, "Danielli_2022/MuSClike.xlsx"), col_names=FALSE)
Danielli_MuSC <- unlist(Danielli_MuSC)
Danielli_cycling <- read_excel(file.path(genelist_dir, "Danielli_2022/Cycling.xlsx"), col_names=FALSE)
Danielli_cycling <- unlist(Danielli_cycling)
Danielli_differentiated <- read_excel(file.path(genelist_dir, "Danielli_2022/Differentiated.xlsx"), col_names=FALSE)
Danielli_differentiated <- unlist(Danielli_differentiated)

Patel_mesoderm <- read_excel(file.path(genelist_dir, "Patel_2022/Mesoderm.xlsx"), col_names=FALSE)
Patel_mesoderm <- unlist(Patel_mesoderm)
Patel_myoblast <- read_excel(file.path(genelist_dir, "Patel_2022/Myoblast.xlsx"), col_names=FALSE)
Patel_myoblast <- unlist(Patel_myoblast)
Patel_myocyte <- read_excel(file.path(genelist_dir, "Patel_2022/Myocyte.xlsx"), col_names=FALSE)
Patel_myocyte <- unlist(Patel_myocyte)


## RMS atlas scores (identified in this new analysis)
Common_differentiated <- read_excel(file.path(genelist_dir, "Common_differentiated.xlsx"), col_names=FALSE)
Common_differentiated <- unlist(Common_differentiated)
Common_proliferative <- read_excel(file.path(genelist_dir,"Common_proliferative.xlsx"), col_names=FALSE)
Common_proliferative <- unlist(Common_proliferative)
Common_Stemcell <- read_excel(file.path(genelist_dir, "Common_Stemcell.xlsx"), col_names=FALSE)
Common_Stemcell <- unlist(Common_Stemcell)



lst <- list(Common_Stemcell, Wei_mesenchymal, Danielli_MuSC, Patel_mesoderm,
            Common_proliferative, Wei_proliferative, Danielli_cycling, Patel_myoblast,
            Common_differentiated,  Wei_muscle, Danielli_differentiated, Patel_myocyte)


#names(lst) <- c('Wei_mesenchymal', 'Wei_proliferative', 'Wei_muscle',
               'Danielli_MuSC',  'Danielli_cycling', 'Danielli_differentiated',
               'Patel_mesoderm', 'Patel_myoblast', 'Patel_myocyte',
               'Common_differentiated', 'Common_proliferative', 'Common_Stemcell')

# Function to calculate Jaccard distance
siml  <- function(x,y) {
  1-length(intersect(lst[[x]],lst[[y]]))/length(union(lst[[x]],lst[[y]]))
}

siml  <- function(x,y) {
  length(intersect(lst[[x]],lst[[y]]))/length(union(lst[[x]],lst[[y]]))
}

# get a data frame of all factor combinations
z <- expand.grid(x=1:length(lst), y=1:length(lst))

# calculate Jaccard distances between each pair of gene sets
result <- mapply(siml,z$x,z$y)
# Turn into matrix
dim(result) <- c(length(lst),length(lst))
# Convert to distance
dists <- as.dist(result)
# Cluster a plot a dendrogram
plot(hclust(dists))

pheatmap(dists)
pheatmap(result)

jaccard <- function(x, y) {
  intersection = length(intersect(lst[[x]], lst[[y]]))
  union = length(lst[[x]]) + length(lst[[y]]) - intersection
  return (intersection/union)
}

result <- mapply(jaccard,z$x,z$y)


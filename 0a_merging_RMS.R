rm(list = ls())

library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(readxl)
library(writexl)
library(data.table)
library(tidyverse)
library(viridis)
library(future)
library(paletteer)
plan("multiprocess", workers = 256)
options(future.globals.maxSize = 256000 * 1024^2)

##################################################################
############ (1) Load datasets Wei et al. & clean ##############
##################################################################

# Load the RMS datasets Langenau

      # 20082
      Langenau_20082 <- readRDS("/mnt/Sara/RMS/Langenau/20082_hg19/20082.rds")
      DefaultAssay(Langenau_20082) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_20082 <- DietSeurat(Langenau_20082, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_20082 <- Langenau_20082[, sample(colnames(Langenau_20082), size = 1500, replace=F)]
      rm(Langenau_20082)
      
      # 20696
      Langenau_20696 <- readRDS("/mnt/Sara/RMS/Langenau/20696_hg19/20696.rds")
      DefaultAssay(Langenau_20696) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_20696 <- DietSeurat(Langenau_20696, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_20696 <- Langenau_20696[, sample(colnames(Langenau_20696), size = 1676, replace=F)]
      rm(Langenau_20696)
      
      # 21202
      Langenau_21202 <- readRDS("/mnt/Sara/RMS/Langenau/21202_hg19/21202.rds")
      DefaultAssay(Langenau_21202) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_21202 <- DietSeurat(Langenau_21202, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_21202 <- Langenau_21202[, sample(colnames(Langenau_21202), size = 1500, replace=F)]
      rm(Langenau_21202)
      
      # 29806
      Langenau_29806 <- readRDS("/mnt/Sara/RMS/Langenau/29806_hg19/29806.rds")
      DefaultAssay(Langenau_29806) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_29806 <- DietSeurat(Langenau_29806, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_29806 <- Langenau_29806[, sample(colnames(Langenau_29806), size = 1500, replace=F)]
      rm(Langenau_29806)
      
      # Mast39
      Langenau_Mast39 <- readRDS("/mnt/Sara/RMS/Langenau/MAST39_hg19_mm10/Mast39.rds")
      DefaultAssay(Langenau_Mast39) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast39 <- DietSeurat(Langenau_Mast39, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast39 <- Langenau_Mast39[, sample(colnames(Langenau_Mast39), size = 1500, replace=F)]
      rm(Langenau_Mast39)
      
      # Mast85 r1
      Langenau_Mast85_r1 <- readRDS("/mnt/Sara/RMS/Langenau/MAST85-r1_hg19_mm10/Mast85r1.rds")
      DefaultAssay(Langenau_Mast85_r1) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast85_r1 <- DietSeurat(Langenau_Mast85_r1, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast85_r1 <- Langenau_Mast85_r1[, sample(colnames(Langenau_Mast85_r1), size = 1500, replace=F)]
      rm(Langenau_Mast85_r1)
      
      
      # Mast85 r2
      Langenau_Mast85_r2 <- readRDS("/mnt/Sara/RMS/Langenau/MAST85-r2_hg19_mm10/Mast85r2.rds")
      DefaultAssay(Langenau_Mast85_r2) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast85_r2 <- DietSeurat(Langenau_Mast85_r2, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast85_r2 <- Langenau_Mast85_r2[, sample(colnames(Langenau_Mast85_r2), size = 1500, replace=F)]
      rm(Langenau_Mast85_r2)
      
      # Mast85 r2-SC
      Langenau_Mast85_r2_SC <- readRDS("/mnt/Sara/RMS/Langenau/MAST85-r2-SC_hg19_mm10/Mast85r2_SC.rds")
      DefaultAssay(Langenau_Mast85_r2_SC) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast85_r2_SC <- DietSeurat(Langenau_Mast85_r2_SC, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast85_r2_SC <- Langenau_Mast85_r2_SC[, sample(colnames(Langenau_Mast85_r2_SC), size = 1500, replace=F)]
      rm(Langenau_Mast85_r2_SC)
      
      # Mast95
      Langenau_Mast95 <- readRDS("/mnt/Sara/RMS/Langenau/MAST95_hg19_mm10/Mast95.rds")
      DefaultAssay(Langenau_Mast95) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast95 <- DietSeurat(Langenau_Mast95, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast95 <- Langenau_Mast95[, sample(colnames(Langenau_Mast95), size = 1500, replace=F)]
      rm(Langenau_Mast95)
      
      # Mast111
      Langenau_Mast111 <- readRDS("/mnt/Sara/RMS/Langenau/MAST111_hg19_mm10/Mast111.rds")
      DefaultAssay(Langenau_Mast111) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast111 <- DietSeurat(Langenau_Mast111, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast111 <- Langenau_Mast111[, sample(colnames(Langenau_Mast111), size = 1500, replace=F)]
      rm(Langenau_Mast111)
      
      # Mast118
      Langenau_Mast118 <- readRDS("/mnt/Sara/RMS/Langenau/MAST118_hg19_mm10/Mast118.rds")
      DefaultAssay(Langenau_Mast118) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast118 <- DietSeurat(Langenau_Mast118, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast118 <- Langenau_Mast118[, sample(colnames(Langenau_Mast118), size = 1500, replace=F)]
      rm(Langenau_Mast118)
      
      # Mast139
      Langenau_Mast139 <- readRDS("/mnt/Sara/RMS/Langenau/MAST139_hg19_mm10/Mast139.rds")
      DefaultAssay(Langenau_Mast139) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast139 <- DietSeurat(Langenau_Mast139, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast139 <- Langenau_Mast139[, sample(colnames(Langenau_Mast139), size = 1500, replace=F)]
      rm(Langenau_Mast139)
      
      # Mast139-SC
      Langenau_Mast139_SC <- readRDS("/mnt/Sara/RMS/Langenau/MAST139-SC_hg19_mm10/Mast139_SC.rds")
      DefaultAssay(Langenau_Mast139_SC) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_Mast139_SC <- DietSeurat(Langenau_Mast139_SC, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_Mast139_SC <- Langenau_Mast139_SC[, sample(colnames(Langenau_Mast139_SC), size = 1500, replace=F)]
      rm(Langenau_Mast139_SC)
      
      
      # MSK72117
      Langenau_MSK72117 <- readRDS("/mnt/Sara/RMS/Langenau/MSK72117_hg19_mm10/MSK72117.rds")
      DefaultAssay(Langenau_MSK72117) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_MSK72117 <- DietSeurat(Langenau_MSK72117, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_MSK72117 <- Langenau_MSK72117[, sample(colnames(Langenau_MSK72117), size = 1500, replace=F)]
      rm(Langenau_MSK72117)
      
      # MSK72117-SC
      Langenau_MSK72117_SC <- readRDS("/mnt/Sara/RMS/Langenau/MSK72117-SC_hg19_mm10/C12SC2_seurat-object.rds")
      DefaultAssay(Langenau_MSK72117_SC) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_MSK72117_SC <- DietSeurat(Langenau_MSK72117_SC, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_MSK72117_SC <- Langenau_MSK72117_SC[, sample(colnames(Langenau_MSK72117_SC), size = 1500, replace=F)]
      rm(Langenau_MSK72117_SC)
      
      
      # MSK74711
      Langenau_MSK74711 <- readRDS("/mnt/Sara/RMS/Langenau/MSK74711_hg19_mm10/MSK74711.rds")
      DefaultAssay(Langenau_MSK74711) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_MSK74711 <- DietSeurat(Langenau_MSK74711, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_MSK74711 <- Langenau_MSK74711[, sample(colnames(Langenau_MSK74711), size = 1500, replace=F)]
      rm(Langenau_MSK74711)
      
      # MSK82489
      Langenau_MSK82489 <- readRDS("/mnt/Sara/RMS/Langenau/MSK82489_hg19_mm10/MSK82489.rds")
      DefaultAssay(Langenau_MSK82489) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_MSK82489 <- DietSeurat(Langenau_MSK82489, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_MSK82489 <- Langenau_MSK82489[, sample(colnames(Langenau_MSK82489), size = 1385, replace=F)]
      rm(Langenau_MSK82489)
      
      # RD
      Langenau_RD <- readRDS("/mnt/Sara/RMS/Langenau/RD_hg19/RD.rds")
      DefaultAssay(Langenau_RD) <- "RNA"
      ## Reduce dimensions of objects
      Langenau_RD <- DietSeurat(Langenau_RD, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      # Subset the dataset to 500 cells
      set.seed(111)
      sm_Langenau_RD <- Langenau_RD[, sample(colnames(Langenau_RD), size = 1500, replace=F)]
      rm(Langenau_RD)
      


##################################################################
############ (2) Load datasets Patel et al. & clean ##############
##################################################################
      
      
      
    # SJRHB011_X
    SJRHB011_X_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB011_X_sc.Rds")
    DefaultAssay(SJRHB011_X_scRNAseq) <- "RNA"
    
    ## Reduce dimensions of objects
    SJRHB011_X_scRNAseq <- DietSeurat(SJRHB011_X_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
   
     ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
    SJRHB011_X_scRNAseq <- SJRHB011_X_scRNAseq[rownames(SJRHB011_X_scRNAseq) %like% "hg19-", ]

    # select only tumor cells from Patel dataset
    SJRHB011_X_scRNAseq <- subset(SJRHB011_X_scRNAseq, cell.cluster.ids == 'Tumor')
    
    ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB011_X_scRNAseq.data <- GetAssayData(object = SJRHB011_X_scRNAseq[["RNA"]], slot = "counts")
      SJRHB011_X_scRNAseq.data <- as.matrix(SJRHB011_X_scRNAseq.data)
    
      # Rename rows
      rownames(SJRHB011_X_scRNAseq.data) <- str_remove_all(rownames(SJRHB011_X_scRNAseq), "hg19-")
    
      # Generate new Seurat object.
      new.SJRHB011_X_scRNAseq <- CreateSeuratObject(
        SJRHB011_X_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB011_X_scRNAseq <- new.SJRHB011_X_scRNAseq
 
      # Subset the dataset to 500 cells
      set.seed(111)
      SJRHB011_X_scRNAseq_small <- SJRHB011_X_scRNAseq[, sample(colnames(SJRHB011_X_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB011_X_scRNAseq)
      
   # SJRHB012_Y  
    SJRHB012_Y_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB012_Y_sc.Rds")
    DefaultAssay(SJRHB012_Y_scRNAseq) <- "RNA"
      
    ## Reduce dimensions of objects
    SJRHB012_Y_scRNAseq <- DietSeurat(SJRHB012_Y_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
    ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
    SJRHB012_Y_scRNAseq <- SJRHB012_Y_scRNAseq[rownames(SJRHB012_Y_scRNAseq) %like% "hg19-", ]
      
    # select only tumor cells from Patel dataset
    SJRHB012_Y_scRNAseq <- subset(SJRHB012_Y_scRNAseq, cell.cluster.ids == 'Tumor')
      
    ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB012_Y_scRNAseq.data <- GetAssayData(object = SJRHB012_Y_scRNAseq[["RNA"]], slot = "counts")
      SJRHB012_Y_scRNAseq.data <- as.matrix(SJRHB012_Y_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB012_Y_scRNAseq.data) <- str_remove_all(rownames(SJRHB012_Y_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB012_Y_scRNAseq <- CreateSeuratObject(
        SJRHB012_Y_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB012_Y_scRNAseq <- new.SJRHB012_Y_scRNAseq  
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB012_Y_scRNAseq_small <- SJRHB012_Y_scRNAseq[, sample(colnames(SJRHB012_Y_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB012_Y_scRNAseq)
      
   # SJRHB012_Z  
    SJRHB012_Z_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB012_Z_sc.Rds")
    DefaultAssay(SJRHB012_Z_scRNAseq) <- "RNA"
      
    ## Reduce dimensions of objects
    SJRHB012_Z_scRNAseq <-  DietSeurat(SJRHB012_Z_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
    ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
    SJRHB012_Z_scRNAseq <- SJRHB012_Z_scRNAseq[rownames(SJRHB012_Z_scRNAseq) %like% "hg19-", ]
      
    # select only tumor cells from Patel dataset
    SJRHB012_Z_scRNAseq <- subset(SJRHB012_Z_scRNAseq, cell.cluster.ids == 'Tumor')
      
    ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB012_Z_scRNAseq.data <- GetAssayData(object = SJRHB012_Z_scRNAseq[["RNA"]], slot = "counts")
      SJRHB012_Z_scRNAseq.data <- as.matrix(SJRHB012_Z_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB012_Z_scRNAseq.data) <- str_remove_all(rownames(SJRHB012_Z_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB012_Z_scRNAseq <- CreateSeuratObject(
        SJRHB012_Z_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB012_Z_scRNAseq <- new.SJRHB012_Z_scRNAseq  
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB012_Z_scRNAseq_small <- SJRHB012_Z_scRNAseq[, sample(colnames(SJRHB012_Z_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB012_Z_scRNAseq)
      
   
  # SJRHB000026_X1_scRNAseq  
    SJRHB000026_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB000026_X1_sc.Rds")
    DefaultAssay(SJRHB000026_X1_scRNAseq) <- "RNA"
      
    ## Reduce dimensions of objects
    SJRHB000026_X1_scRNAseq <- DietSeurat(SJRHB000026_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB000026_X1_scRNAseq <- SJRHB000026_X1_scRNAseq[rownames(SJRHB000026_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB000026_X1_scRNAseq <- subset(SJRHB000026_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB000026_X1_scRNAseq.data <- GetAssayData(object = SJRHB000026_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB000026_X1_scRNAseq.data <- as.matrix(SJRHB000026_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB000026_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB000026_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB000026_X1_scRNAseq <- CreateSeuratObject(
        SJRHB000026_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB000026_X1_scRNAseq <- new.SJRHB000026_X1_scRNAseq     
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB000026_X1_scRNAseq_small <- SJRHB000026_X1_scRNAseq[, sample(colnames(SJRHB000026_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB000026_X1_scRNAseq)
      
 
# SJRHB000026_R2
      SJRHB000026_R2_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB000026_R2_sn.Rds")
      DefaultAssay(SJRHB000026_R2_snRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB000026_R2_snRNAseq <- DietSeurat(SJRHB000026_R2_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB000026_R2_snRNAseq <- SJRHB000026_R2_snRNAseq[rownames(SJRHB000026_R2_snRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB000026_R2_snRNAseq <- subset(SJRHB000026_R2_snRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB000026_R2_snRNAseq.data <- GetAssayData(object = SJRHB000026_R2_snRNAseq[["RNA"]], slot = "counts")
      SJRHB000026_R2_snRNAseq.data <- as.matrix(SJRHB000026_R2_snRNAseq.data)
      
      # Rename rows
      rownames(SJRHB000026_R2_snRNAseq.data) <- str_remove_all(rownames(SJRHB000026_R2_snRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB000026_R2_snRNAseq <- CreateSeuratObject(
        SJRHB000026_R2_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB000026_R2_snRNAseq <- new.SJRHB000026_R2_snRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB000026_R2_snRNAseq_small <- SJRHB000026_R2_snRNAseq[, sample(colnames(SJRHB000026_R2_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB000026_R2_snRNAseq)
 
      
# SJRHB000026_R3
      SJRHB000026_R3_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB000026_R3_sn.Rds")
      DefaultAssay(SJRHB000026_R3_snRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB000026_R3_snRNAseq <- DietSeurat(SJRHB000026_R3_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB000026_R3_snRNAseq <- SJRHB000026_R3_snRNAseq[rownames(SJRHB000026_R3_snRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB000026_R3_snRNAseq <- subset(SJRHB000026_R3_snRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB000026_R3_snRNAseq.data <- GetAssayData(object = SJRHB000026_R3_snRNAseq[["RNA"]], slot = "counts")
      SJRHB000026_R3_snRNAseq.data <- as.matrix(SJRHB000026_R3_snRNAseq.data)
      
      # Rename rows
      rownames(SJRHB000026_R3_snRNAseq.data) <- str_remove_all(rownames(SJRHB000026_R3_snRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB000026_R3_snRNAseq <- CreateSeuratObject(
        SJRHB000026_R3_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB000026_R3_snRNAseq <- new.SJRHB000026_R3_snRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB000026_R3_snRNAseq_small <- SJRHB000026_R3_snRNAseq[, sample(colnames(SJRHB000026_R3_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB000026_R3_snRNAseq)
      

      
# SJRHB000026_X2_scRNAseq 
      SJRHB000026_X2_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB000026_X2_sc.Rds")
      DefaultAssay(SJRHB000026_X2_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB000026_X2_scRNAseq <- DietSeurat(SJRHB000026_X2_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB000026_X2_scRNAseq <- SJRHB000026_X2_scRNAseq[rownames(SJRHB000026_X2_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB000026_X2_scRNAseq <- subset(SJRHB000026_X2_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB000026_X2_scRNAseq.data <- GetAssayData(object = SJRHB000026_X2_scRNAseq[["RNA"]], slot = "counts")
      SJRHB000026_X2_scRNAseq.data <- as.matrix(SJRHB000026_X2_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB000026_X2_scRNAseq.data) <- str_remove_all(rownames(SJRHB000026_X2_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB000026_X2_scRNAseq <- CreateSeuratObject(
        SJRHB000026_X2_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB000026_X2_scRNAseq <- new.SJRHB000026_X2_scRNAseq     
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB000026_X2_scRNAseq_small <- SJRHB000026_X2_scRNAseq[, sample(colnames(SJRHB000026_X2_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB000026_X2_scRNAseq)
      
      
      
# SJRHB000026_X1_scRNAseq  
      SJRHB000026_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB000026_X1_sc.Rds")
      DefaultAssay(SJRHB000026_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB000026_X1_scRNAseq <- DietSeurat(SJRHB000026_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB000026_X1_scRNAseq <- SJRHB000026_X1_scRNAseq[rownames(SJRHB000026_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB000026_X1_scRNAseq <- subset(SJRHB000026_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB000026_X1_scRNAseq.data <- GetAssayData(object = SJRHB000026_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB000026_X1_scRNAseq.data <- as.matrix(SJRHB000026_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB000026_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB000026_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB000026_X1_scRNAseq <- CreateSeuratObject(
        SJRHB000026_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB000026_X1_scRNAseq <- new.SJRHB000026_X1_scRNAseq     
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB000026_X1_scRNAseq_small <- SJRHB000026_X1_scRNAseq[, sample(colnames(SJRHB000026_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB000026_X1_scRNAseq)
      
      
   # SJRHB010468_X1_scRNAseq  
      SJRHB010468_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB010468_X1_sc.Rds")
      DefaultAssay(SJRHB010468_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB010468_X1_scRNAseq <- DietSeurat(SJRHB010468_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB010468_X1_scRNAseq <- SJRHB010468_X1_scRNAseq[rownames(SJRHB010468_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB010468_X1_scRNAseq <- subset(SJRHB010468_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB010468_X1_scRNAseq.data <- GetAssayData(object = SJRHB010468_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB010468_X1_scRNAseq.data <- as.matrix(SJRHB010468_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB010468_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB010468_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB010468_X1_scRNAseq <- CreateSeuratObject(
        SJRHB010468_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB010468_X1_scRNAseq <- new.SJRHB010468_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB010468_X1_scRNAseq_small <- SJRHB010468_X1_scRNAseq[, sample(colnames(SJRHB010468_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB010468_X1_scRNAseq)
      
    
  #SJRHB010927_D1_snRNAseq  
      SJRHB010927_D1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB010927_D1_sn.Rds")
      DefaultAssay(SJRHB010927_D1_snRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB010927_D1_snRNAseq <- DietSeurat(SJRHB010927_D1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB010927_D1_snRNAseq <- SJRHB010927_D1_snRNAseq[rownames(SJRHB010927_D1_snRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB010927_D1_snRNAseq <- subset(SJRHB010927_D1_snRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB010927_D1_snRNAseq.data <- GetAssayData(object = SJRHB010927_D1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB010927_D1_snRNAseq.data <- as.matrix(SJRHB010927_D1_snRNAseq.data)
      
      # Rename rows
      rownames(SJRHB010927_D1_snRNAseq.data) <- str_remove_all(rownames(SJRHB010927_D1_snRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB010927_D1_snRNAseq <- CreateSeuratObject(
        SJRHB010927_D1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB010927_D1_snRNAseq <- new.SJRHB010927_D1_snRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB010927_D1_snRNAseq_small <- SJRHB010927_D1_snRNAseq[, sample(colnames(SJRHB010927_D1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB010927_D1_snRNAseq)
   
   
#SJRHB010927_X1_scRNAseq  
      SJRHB010927_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB010927_X1_sc.Rds")
      DefaultAssay(SJRHB010927_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB010927_X1_scRNAseq <- DietSeurat(SJRHB010927_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB010927_X1_scRNAseq <- SJRHB010927_X1_scRNAseq[rownames(SJRHB010927_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB010927_X1_scRNAseq <- subset(SJRHB010927_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB010927_X1_scRNAseq.data <- GetAssayData(object = SJRHB010927_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB010927_X1_scRNAseq.data <- as.matrix(SJRHB010927_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB010927_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB010927_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB010927_X1_scRNAseq <- CreateSeuratObject(
        SJRHB010927_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB010927_X1_scRNAseq <- new.SJRHB010927_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB010927_X1_scRNAseq_small <- SJRHB010927_X1_scRNAseq[, sample(colnames(SJRHB010927_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB010927_X1_scRNAseq)
      
      
      
         
   # SJRHB010928_X1_scRNAseq  
      SJRHB010928_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB010928_X1_sc.Rds")
      DefaultAssay(SJRHB010928_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB010928_X1_scRNAseq <- DietSeurat(SJRHB010928_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB010928_X1_scRNAseq <- SJRHB010928_X1_scRNAseq[rownames(SJRHB010928_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB010928_X1_scRNAseq <- subset(SJRHB010928_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB010928_X1_scRNAseq.data <- GetAssayData(object = SJRHB010928_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB010928_X1_scRNAseq.data <- as.matrix(SJRHB010928_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB010928_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB010928_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB010928_X1_scRNAseq <- CreateSeuratObject(
        SJRHB010928_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB010928_X1_scRNAseq <- new.SJRHB010928_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB010928_X1_scRNAseq_small <- SJRHB010928_X1_scRNAseq[, sample(colnames(SJRHB010928_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB010928_X1_scRNAseq)   

      
# SJRHB010928_R1 
      SJRHB010928_R1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB010928_R1_sn.Rds")
      DefaultAssay(SJRHB010928_R1_snRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB010928_R1_snRNAseq <- DietSeurat(SJRHB010928_R1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB010928_R1_snRNAseq <- SJRHB010928_R1_snRNAseq[rownames(SJRHB010928_R1_snRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB010928_R1_snRNAseq <- subset(SJRHB010928_R1_snRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB010928_R1_snRNAseq.data <- GetAssayData(object = SJRHB010928_R1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB010928_R1_snRNAseq.data <- as.matrix(SJRHB010928_R1_snRNAseq.data)
      
      # Rename rows
      rownames(SJRHB010928_R1_snRNAseq.data) <- str_remove_all(rownames(SJRHB010928_R1_snRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB010928_R1_snRNAseq <- CreateSeuratObject(
        SJRHB010928_R1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB010928_R1_snRNAseq <- new.SJRHB010928_R1_snRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB010928_R1_snRNAseq_small <- SJRHB010928_R1_snRNAseq[, sample(colnames(SJRHB010928_R1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB010928_R1_snRNAseq)  
      
      
            

  
   # SJRHB012405_X1_scRNAseq  
      SJRHB012405_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB012405_X1_sc.Rds")
      DefaultAssay(SJRHB012405_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB012405_X1_scRNAseq <-  DietSeurat(SJRHB012405_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB012405_X1_scRNAseq <- SJRHB012405_X1_scRNAseq[rownames(SJRHB012405_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB012405_X1_scRNAseq <- subset(SJRHB012405_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB012405_X1_scRNAseq.data <- GetAssayData(object = SJRHB012405_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB012405_X1_scRNAseq.data <- as.matrix(SJRHB012405_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB012405_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB012405_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB012405_X1_scRNAseq <- CreateSeuratObject(
        SJRHB012405_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB012405_X1_scRNAseq <- new.SJRHB012405_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB012405_X1_scRNAseq_small <- SJRHB012405_X1_scRNAseq[, sample(colnames(SJRHB012405_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB012405_X1_scRNAseq)   
      
      
      
      
 # SJRHB012405_D1
      SJRHB012405_D1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB012405_D1_sn.Rds")
      DefaultAssay(SJRHB012405_D1_snRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB012405_D1_snRNAseq <-  DietSeurat(SJRHB012405_D1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB012405_D1_snRNAseq <- SJRHB012405_D1_snRNAseq[rownames(SJRHB012405_D1_snRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB012405_D1_snRNAseq <- subset(SJRHB012405_D1_snRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB012405_D1_snRNAseq.data <- GetAssayData(object = SJRHB012405_D1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB012405_D1_snRNAseq.data <- as.matrix(SJRHB012405_D1_snRNAseq.data)
      
      # Rename rows
      rownames(SJRHB012405_D1_snRNAseq.data) <- str_remove_all(rownames(SJRHB012405_D1_snRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB012405_D1_snRNAseq <- CreateSeuratObject(
        SJRHB012405_D1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB012405_D1_snRNAseq <- new.SJRHB012405_D1_snRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB012405_D1_snRNAseq_small <- SJRHB012405_D1_snRNAseq[, sample(colnames(SJRHB012405_D1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB012405_D1_snRNAseq)   
      
      
      
   #SJRHB013757_X1_scRNAseq  
      SJRHB013757_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB013757_X1_sc.Rds")
      DefaultAssay(SJRHB013757_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB013757_X1_scRNAseq <- DietSeurat(SJRHB013757_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013757_X1_scRNAseq <- SJRHB013757_X1_scRNAseq[rownames(SJRHB013757_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB013757_X1_scRNAseq <- subset(SJRHB013757_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013757_X1_scRNAseq.data <- GetAssayData(object = SJRHB013757_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB013757_X1_scRNAseq.data <- as.matrix(SJRHB013757_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB013757_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB013757_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB013757_X1_scRNAseq <- CreateSeuratObject(
        SJRHB013757_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB013757_X1_scRNAseq <- new.SJRHB013757_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013757_X1_scRNAseq_small <- SJRHB013757_X1_scRNAseq[, sample(colnames(SJRHB013757_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB013757_X1_scRNAseq) 
      
   #SJRHB013758_X1_scRNAseq 
      SJRHB013758_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB013758_X1_sc.Rds")
      DefaultAssay(SJRHB013758_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB013758_X1_scRNAseq <- DietSeurat(SJRHB013758_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013758_X1_scRNAseq <- SJRHB013758_X1_scRNAseq[rownames(SJRHB013758_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB013758_X1_scRNAseq <- subset(SJRHB013758_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013758_X1_scRNAseq.data <- GetAssayData(object = SJRHB013758_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB013758_X1_scRNAseq.data <- as.matrix(SJRHB013758_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB013758_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB013758_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB013758_X1_scRNAseq <- CreateSeuratObject(
        SJRHB013758_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB013758_X1_scRNAseq <- new.SJRHB013758_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013758_X1_scRNAseq_small <- SJRHB013758_X1_scRNAseq[, sample(colnames(SJRHB013758_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB013758_X1_scRNAseq)      
      
   #SJRHB013758_X2_scRNAseq 
      SJRHB013758_X2_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB013758_X2_sc.Rds")
      DefaultAssay(SJRHB013758_X2_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      SJRHB013758_X2_scRNAseq <- DietSeurat(SJRHB013758_X2_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013758_X2_scRNAseq <- SJRHB013758_X2_scRNAseq[rownames(SJRHB013758_X2_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB013758_X2_scRNAseq <- subset(SJRHB013758_X2_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013758_X2_scRNAseq.data <- GetAssayData(object = SJRHB013758_X2_scRNAseq[["RNA"]], slot = "counts")
      SJRHB013758_X2_scRNAseq.data <- as.matrix(SJRHB013758_X2_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB013758_X2_scRNAseq.data) <- str_remove_all(rownames(SJRHB013758_X2_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB013758_X2_scRNAseq <- CreateSeuratObject(
        SJRHB013758_X2_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB013758_X2_scRNAseq <- new.SJRHB013758_X2_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013758_X2_scRNAseq_small <- SJRHB013758_X2_scRNAseq[, sample(colnames(SJRHB013758_X2_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB013758_X2_scRNAseq)      
    
   #SJRHB013759_X14_scRNAseq 
      SJRHB013759_X14_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB013759_X14_sc.Rds")
      DefaultAssay(SJRHB013759_X14_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB013759_X14_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013759_X14_scRNAseq <- SJRHB013759_X14_scRNAseq[rownames(SJRHB013759_X14_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB013759_X14_scRNAseq <- subset(SJRHB013759_X14_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013759_X14_scRNAseq.data <- GetAssayData(object = SJRHB013759_X14_scRNAseq[["RNA"]], slot = "counts")
      SJRHB013759_X14_scRNAseq.data <- as.matrix(SJRHB013759_X14_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB013759_X14_scRNAseq.data) <- str_remove_all(rownames(SJRHB013759_X14_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB013759_X14_scRNAseq <- CreateSeuratObject(
        SJRHB013759_X14_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB013759_X14_scRNAseq <- new.SJRHB013759_X14_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013759_X14_scRNAseq_small <- SJRHB013759_X14_scRNAseq[, sample(colnames(SJRHB013759_X14_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB013759_X14_scRNAseq) 
      
  #SJRHB013759_X15_scRNAseq 
      SJRHB013759_X15_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB013759_X15_sc.Rds")
      DefaultAssay(SJRHB013759_X15_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB013759_X15_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013759_X15_scRNAseq <- SJRHB013759_X15_scRNAseq[rownames(SJRHB013759_X15_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB013759_X15_scRNAseq <- subset(SJRHB013759_X15_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013759_X15_scRNAseq.data <- GetAssayData(object = SJRHB013759_X15_scRNAseq[["RNA"]], slot = "counts")
      SJRHB013759_X15_scRNAseq.data <- as.matrix(SJRHB013759_X15_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB013759_X15_scRNAseq.data) <- str_remove_all(rownames(SJRHB013759_X15_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB013759_X15_scRNAseq <- CreateSeuratObject(
        SJRHB013759_X15_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB013759_X15_scRNAseq <- new.SJRHB013759_X15_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013759_X15_scRNAseq_small <- SJRHB013759_X15_scRNAseq[, sample(colnames(SJRHB013759_X15_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB013759_X15_scRNAseq) 
    
  #SJRHB030680_X1_scRNAseq 
      SJRHB030680_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB030680_X1_sc.Rds")
      DefaultAssay(SJRHB030680_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB030680_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB030680_X1_scRNAseq <- SJRHB030680_X1_scRNAseq[rownames(SJRHB030680_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB030680_X1_scRNAseq <- subset(SJRHB030680_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB030680_X1_scRNAseq.data <- GetAssayData(object = SJRHB030680_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB030680_X1_scRNAseq.data <- as.matrix(SJRHB030680_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB030680_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB030680_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB030680_X1_scRNAseq <- CreateSeuratObject(
        SJRHB030680_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB030680_X1_scRNAseq <- new.SJRHB030680_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB030680_X1_scRNAseq_small <- SJRHB030680_X1_scRNAseq[, sample(colnames(SJRHB030680_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB030680_X1_scRNAseq) 
      

      
  #SJRHB031320_X1_scRNAseq 
      SJRHB031320_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB031320_X1_sc.Rds")
      DefaultAssay(SJRHB031320_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB031320_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB031320_X1_scRNAseq <- SJRHB031320_X1_scRNAseq[rownames(SJRHB031320_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB031320_X1_scRNAseq <- subset(SJRHB031320_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB031320_X1_scRNAseq.data <- GetAssayData(object = SJRHB031320_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB031320_X1_scRNAseq.data <- as.matrix(SJRHB031320_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB031320_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB031320_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB031320_X1_scRNAseq <- CreateSeuratObject(
        SJRHB031320_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB031320_X1_scRNAseq <- new.SJRHB031320_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB031320_X1_scRNAseq_small <- SJRHB031320_X1_scRNAseq[, sample(colnames(SJRHB031320_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB031320_X1_scRNAseq) 
      
  #SJRHB046156_X1_scRNAseq 
      SJRHB046156_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB046156_X1_sc.Rds")
      DefaultAssay(SJRHB046156_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB046156_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB046156_X1_scRNAseq <- SJRHB046156_X1_scRNAseq[rownames(SJRHB046156_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB046156_X1_scRNAseq <- subset(SJRHB046156_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB046156_X1_scRNAseq.data <- GetAssayData(object = SJRHB046156_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB046156_X1_scRNAseq.data <- as.matrix(SJRHB046156_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB046156_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB046156_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB046156_X1_scRNAseq <- CreateSeuratObject(
        SJRHB046156_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB046156_X1_scRNAseq <- new.SJRHB046156_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB046156_X1_scRNAseq_small <- SJRHB046156_X1_scRNAseq[, sample(colnames(SJRHB046156_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB046156_X1_scRNAseq)       
    
  #SJRHB049189_X1_scRNAseq 
      SJRHB049189_X1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB049189_X1_sc.Rds")
      DefaultAssay(SJRHB049189_X1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB049189_X1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB049189_X1_scRNAseq <- SJRHB049189_X1_scRNAseq[rownames(SJRHB049189_X1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB049189_X1_scRNAseq <- subset(SJRHB049189_X1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB049189_X1_scRNAseq.data <- GetAssayData(object = SJRHB049189_X1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB049189_X1_scRNAseq.data <- as.matrix(SJRHB049189_X1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB049189_X1_scRNAseq.data) <- str_remove_all(rownames(SJRHB049189_X1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB049189_X1_scRNAseq <- CreateSeuratObject(
        SJRHB049189_X1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB049189_X1_scRNAseq <- new.SJRHB049189_X1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB049189_X1_scRNAseq_small <- SJRHB049189_X1_scRNAseq[, sample(colnames(SJRHB049189_X1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB049189_X1_scRNAseq) 
 
      
#SJRHB049189_D1
      SJRHB049189_D1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB049189_D1_sn.Rds")
      DefaultAssay(SJRHB049189_D1_snRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB049189_D1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB049189_D1_snRNAseq <- SJRHB049189_D1_snRNAseq[rownames(SJRHB049189_D1_snRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB049189_D1_snRNAseq <- subset(SJRHB049189_D1_snRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB049189_D1_snRNAseq.data <- GetAssayData(object = SJRHB049189_D1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB049189_D1_snRNAseq.data <- as.matrix(SJRHB049189_D1_snRNAseq.data)
      
      # Rename rows
      rownames(SJRHB049189_D1_snRNAseq.data) <- str_remove_all(rownames(SJRHB049189_D1_snRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB049189_D1_snRNAseq <- CreateSeuratObject(
        SJRHB049189_D1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB049189_D1_snRNAseq <- new.SJRHB049189_D1_snRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB049189_D1_snRNAseq_small <- SJRHB049189_D1_snRNAseq[, sample(colnames(SJRHB049189_D1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB049189_D1_snRNAseq) 
      
      
           

  #SJRHB011_D_snRNAseq
      SJRHB011_D_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB011_D_sn.Rds")
      DefaultAssay(SJRHB011_D_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB011_D_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB011_D_snRNAseq <- SJRHB011_D_snRNAseq[rownames(SJRHB011_D_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB011_D_snRNAseq <- subset(SJRHB011_D_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB011_D_snRNAseq.data <- GetAssayData(object = SJRHB011_D_snRNAseq[["RNA"]], slot = "counts")
      SJRHB011_D_snRNAseq.data <- as.matrix(SJRHB011_D_snRNAseq.data)

      # Rename rows
      rownames(SJRHB011_D_snRNAseq.data) <- str_remove_all(rownames(SJRHB011_D_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB011_D_snRNAseq <- CreateSeuratObject(
        SJRHB011_D_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB011_D_snRNAseq <- new.SJRHB011_D_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB011_D_snRNAseq_small <- SJRHB011_D_snRNAseq[, sample(colnames(SJRHB011_D_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB011_D_snRNAseq)


  #SJRHB012_R_snRNAseq
      SJRHB012_R_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB012_R_sn.Rds")
      DefaultAssay(SJRHB012_R_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB012_R_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB012_R_snRNAseq <- SJRHB012_R_snRNAseq[rownames(SJRHB012_R_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB012_R_snRNAseq <- subset(SJRHB012_R_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB012_R_snRNAseq.data <- GetAssayData(object = SJRHB012_R_snRNAseq[["RNA"]], slot = "counts")
      SJRHB012_R_snRNAseq.data <- as.matrix(SJRHB012_R_snRNAseq.data)

      # Rename rows
      rownames(SJRHB012_R_snRNAseq.data) <- str_remove_all(rownames(SJRHB012_R_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB012_R_snRNAseq <- CreateSeuratObject(
        SJRHB012_R_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB012_R_snRNAseq <- new.SJRHB012_R_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB012_R_snRNAseq_small <- SJRHB012_R_snRNAseq[, sample(colnames(SJRHB012_R_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB012_R_snRNAseq)

  #SJRHB012_S_snRNAseq
      SJRHB012_S_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB012_S_sn.Rds")
      DefaultAssay(SJRHB012_S_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB012_S_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB012_S_snRNAseq <- SJRHB012_S_snRNAseq[rownames(SJRHB012_S_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB012_S_snRNAseq <- subset(SJRHB012_S_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB012_S_snRNAseq.data <- GetAssayData(object = SJRHB012_S_snRNAseq[["RNA"]], slot = "counts")
      SJRHB012_S_snRNAseq.data <- as.matrix(SJRHB012_S_snRNAseq.data)

      # Rename rows
      rownames(SJRHB012_S_snRNAseq.data) <- str_remove_all(rownames(SJRHB012_S_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB012_S_snRNAseq <- CreateSeuratObject(
        SJRHB012_S_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB012_S_snRNAseq <- new.SJRHB012_S_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB012_S_snRNAseq_small <- SJRHB012_S_snRNAseq[, sample(colnames(SJRHB012_S_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB012_S_snRNAseq)

  #SJRHB010468_D1_snRNAseq
      SJRHB010468_D1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB010468_D1_sn.Rds")
      DefaultAssay(SJRHB010468_D1_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB010468_D1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB010468_D1_snRNAseq <- SJRHB010468_D1_snRNAseq[rownames(SJRHB010468_D1_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB010468_D1_snRNAseq <- subset(SJRHB010468_D1_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB010468_D1_snRNAseq.data <- GetAssayData(object = SJRHB010468_D1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB010468_D1_snRNAseq.data <- as.matrix(SJRHB010468_D1_snRNAseq.data)

      # Rename rows
      rownames(SJRHB010468_D1_snRNAseq.data) <- str_remove_all(rownames(SJRHB010468_D1_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB010468_D1_snRNAseq <- CreateSeuratObject(
        SJRHB010468_D1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB010468_D1_snRNAseq <- new.SJRHB010468_D1_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB010468_D1_snRNAseq_small <- SJRHB010468_D1_snRNAseq[, sample(colnames(SJRHB010468_D1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB010468_D1_snRNAseq)

  #SJRHB013757_D2_snRNAseq
      SJRHB013757_D2_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB013757_D2_sn.Rds")
      DefaultAssay(SJRHB013757_D2_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB013757_D2_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013757_D2_snRNAseq <- SJRHB013757_D2_snRNAseq[rownames(SJRHB013757_D2_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB013757_D2_snRNAseq <- subset(SJRHB013757_D2_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013757_D2_snRNAseq.data <- GetAssayData(object = SJRHB013757_D2_snRNAseq[["RNA"]], slot = "counts")
      SJRHB013757_D2_snRNAseq.data <- as.matrix(SJRHB013757_D2_snRNAseq.data)

      # Rename rows
      rownames(SJRHB013757_D2_snRNAseq.data) <- str_remove_all(rownames(SJRHB013757_D2_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB013757_D2_snRNAseq <- CreateSeuratObject(
        SJRHB013757_D2_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB013757_D2_snRNAseq <- new.SJRHB013757_D2_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013757_D2_snRNAseq_small <- SJRHB013757_D2_snRNAseq[, sample(colnames(SJRHB013757_D2_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB013757_D2_snRNAseq)

  #SJRHB013758_D1_snRNAseq
      SJRHB013758_D1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB013758_D1_sn.Rds")
      DefaultAssay(SJRHB013758_D1_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB013758_D1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013758_D1_snRNAseq <- SJRHB013758_D1_snRNAseq[rownames(SJRHB013758_D1_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB013758_D1_snRNAseq <- subset(SJRHB013758_D1_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013758_D1_snRNAseq.data <- GetAssayData(object = SJRHB013758_D1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB013758_D1_snRNAseq.data <- as.matrix(SJRHB013758_D1_snRNAseq.data)

      # Rename rows
      rownames(SJRHB013758_D1_snRNAseq.data) <- str_remove_all(rownames(SJRHB013758_D1_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB013758_D1_snRNAseq <- CreateSeuratObject(
        SJRHB013758_D1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB013758_D1_snRNAseq <- new.SJRHB013758_D1_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013758_D1_snRNAseq_small <- SJRHB013758_D1_snRNAseq[, sample(colnames(SJRHB013758_D1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB013758_D1_snRNAseq)

  #SJRHB013758_D2_snRNAseq
      SJRHB013758_D2_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB013758_D2_sn.Rds")
      DefaultAssay(SJRHB013758_D2_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB013758_D2_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013758_D2_snRNAseq <- SJRHB013758_D2_snRNAseq[rownames(SJRHB013758_D2_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB013758_D2_snRNAseq <- subset(SJRHB013758_D2_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013758_D2_snRNAseq.data <- GetAssayData(object = SJRHB013758_D2_snRNAseq[["RNA"]], slot = "counts")
      SJRHB013758_D2_snRNAseq.data <- as.matrix(SJRHB013758_D2_snRNAseq.data)

      # Rename rows
      rownames(SJRHB013758_D2_snRNAseq.data) <- str_remove_all(rownames(SJRHB013758_D2_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB013758_D2_snRNAseq <- CreateSeuratObject(
        SJRHB013758_D2_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB013758_D2_snRNAseq <- new.SJRHB013758_D2_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013758_D2_snRNAseq_small <- SJRHB013758_D2_snRNAseq[, sample(colnames(SJRHB013758_D2_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB013758_D2_snRNAseq)

  #SJRHB013759_A1_snRNAseq
      SJRHB013759_A1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB013759_A1_sn.Rds")
      DefaultAssay(SJRHB013759_A1_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB013759_A1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013759_A1_snRNAseq <- SJRHB013759_A1_snRNAseq[rownames(SJRHB013759_A1_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB013759_A1_snRNAseq <- subset(SJRHB013759_A1_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013759_A1_snRNAseq.data <- GetAssayData(object = SJRHB013759_A1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB013759_A1_snRNAseq.data <- as.matrix(SJRHB013759_A1_snRNAseq.data)

      # Rename rows
      rownames(SJRHB013759_A1_snRNAseq.data) <- str_remove_all(rownames(SJRHB013759_A1_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB013759_A1_snRNAseq <- CreateSeuratObject(
        SJRHB013759_A1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB013759_A1_snRNAseq <- new.SJRHB013759_A1_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013759_A1_snRNAseq_small <- SJRHB013759_A1_snRNAseq[, sample(colnames(SJRHB013759_A1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB013759_A1_snRNAseq)


  #SJRHB013759_A2_snRNAseq
      SJRHB013759_A2_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB013759_A2_sn.Rds")
      DefaultAssay(SJRHB013759_A2_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB013759_A2_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB013759_A2_snRNAseq <- SJRHB013759_A2_snRNAseq[rownames(SJRHB013759_A2_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB013759_A2_snRNAseq <- subset(SJRHB013759_A2_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB013759_A2_snRNAseq.data <- GetAssayData(object = SJRHB013759_A2_snRNAseq[["RNA"]], slot = "counts")
      SJRHB013759_A2_snRNAseq.data <- as.matrix(SJRHB013759_A2_snRNAseq.data)

      # Rename rows
      rownames(SJRHB013759_A2_snRNAseq.data) <- str_remove_all(rownames(SJRHB013759_A2_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB013759_A2_snRNAseq <- CreateSeuratObject(
        SJRHB013759_A2_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB013759_A2_snRNAseq <- new.SJRHB013759_A2_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB013759_A2_snRNAseq_small <- SJRHB013759_A2_snRNAseq[, sample(colnames(SJRHB013759_A2_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB013759_A2_snRNAseq)

  #SJRHB046156_A1_snRNAseq
      SJRHB046156_A1_snRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/snRNAseq/SJRHB046156_A1_sn.Rds")
      DefaultAssay(SJRHB046156_A1_snRNAseq) <- "RNA"

      ## Reduce dimensions of objects
      DietSeurat(SJRHB046156_A1_snRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)

      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB046156_A1_snRNAseq <- SJRHB046156_A1_snRNAseq[rownames(SJRHB046156_A1_snRNAseq) %like% "hg19-", ]

      # select only tumor cells from Patel dataset
      SJRHB046156_A1_snRNAseq <- subset(SJRHB046156_A1_snRNAseq, cell.cluster.ids == 'Tumor')

      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB046156_A1_snRNAseq.data <- GetAssayData(object = SJRHB046156_A1_snRNAseq[["RNA"]], slot = "counts")
      SJRHB046156_A1_snRNAseq.data <- as.matrix(SJRHB046156_A1_snRNAseq.data)

      # Rename rows
      rownames(SJRHB046156_A1_snRNAseq.data) <- str_remove_all(rownames(SJRHB046156_A1_snRNAseq), "hg19-")

      # Generate new Seurat object.
      new.SJRHB046156_A1_snRNAseq <- CreateSeuratObject(
        SJRHB046156_A1_snRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )

      SJRHB046156_A1_snRNAseq <- new.SJRHB046156_A1_snRNAseq

      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB046156_A1_snRNAseq_small <- SJRHB046156_A1_snRNAseq[, sample(colnames(SJRHB046156_A1_snRNAseq), size = 1500, replace=F)]
      rm(SJRHB046156_A1_snRNAseq)
      
  #SJRHB031320_D1_scRNAseq
      SJRHB031320_D1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB031320_D1_sc.Rds")
      DefaultAssay(SJRHB031320_D1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB031320_D1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB031320_D1_scRNAseq <- SJRHB031320_D1_scRNAseq[rownames(SJRHB031320_D1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB031320_D1_scRNAseq <- subset(SJRHB031320_D1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB031320_D1_scRNAseq.data <- GetAssayData(object = SJRHB031320_D1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB031320_D1_scRNAseq.data <- as.matrix(SJRHB031320_D1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB031320_D1_scRNAseq.data) <- str_remove_all(rownames(SJRHB031320_D1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB031320_D1_scRNAseq <- CreateSeuratObject(
        SJRHB031320_D1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB031320_D1_scRNAseq <- new.SJRHB031320_D1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB031320_D1_scRNAseq_small <- SJRHB031320_D1_scRNAseq[, sample(colnames(SJRHB031320_D1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB031320_D1_scRNAseq) 

  #SJRHB030680_R1_scRNAseq
      SJRHB030680_R1_scRNAseq <- readRDS("/mnt/Sara/RMS/Patel_2022/SJRHB030680_R1_sc.Rds")
      DefaultAssay(SJRHB030680_R1_scRNAseq) <- "RNA"
      
      ## Reduce dimensions of objects
      DietSeurat(SJRHB030680_R1_scRNAseq, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
      
      ## Remove all rows with mouse genes (mm10-) = select all human genes (hg19-)
      SJRHB030680_R1_scRNAseq <- SJRHB030680_R1_scRNAseq[rownames(SJRHB030680_R1_scRNAseq) %like% "hg19-", ]
      
      # select only tumor cells from Patel dataset
      SJRHB030680_R1_scRNAseq <- subset(SJRHB030680_R1_scRNAseq, cell.cluster.ids == 'Tumor')
      
      ### Remove gene initials "hg19- from Patel"
      # get assay data
      SJRHB030680_R1_scRNAseq.data <- GetAssayData(object = SJRHB030680_R1_scRNAseq[["RNA"]], slot = "counts")
      SJRHB030680_R1_scRNAseq.data <- as.matrix(SJRHB030680_R1_scRNAseq.data)
      
      # Rename rows
      rownames(SJRHB030680_R1_scRNAseq.data) <- str_remove_all(rownames(SJRHB030680_R1_scRNAseq), "hg19-")
      
      # Generate new Seurat object.
      new.SJRHB030680_R1_scRNAseq <- CreateSeuratObject(
        SJRHB030680_R1_scRNAseq.data,
        project = "SeuratProject",
        assay = "RNA",
        min.cells = 0,
        min.features = 0,
        names.field = 1,
        names.delim = "_",
        meta.data = NULL
      )
      
      SJRHB030680_R1_scRNAseq <- new.SJRHB030680_R1_scRNAseq 
      
      # Subset the cells to 1500 cells
      set.seed(111)
      SJRHB030680_R1_scRNAseq_small <- SJRHB030680_R1_scRNAseq[, sample(colnames(SJRHB030680_R1_scRNAseq), size = 1500, replace=F)]
      rm(SJRHB030680_R1_scRNAseq) 
      


##################################################################
############ (3) Load datasets Danielli et al. + Rh41 & clean ####
##################################################################
    
# Load the RMS datasets Danielli et al
    PDX104_TS <- readRDS("/mnt/Sara/RMS/write/PDX104_TS.rds")
    PDX29_TS <- readRDS("/mnt/Sara/RMS/write/PDX29_TS.rds")
    PDX35_TS <- readRDS("/mnt/Sara/RMS/write/PDX35_TS.rds")
    KFR_2 <- readRDS("/mnt/Sara/RMS/write/KFR_2.rds")
    Rh4_2 <- readRDS("/mnt/Sara/RMS/write/Rh4_2.rds")
    RMS_2 <- readRDS("/mnt/Sara/RMS/write/RMS_2.rds")
    Mast118 <- readRDS("/mnt/Sara/RMS/write/Mast118.rds")
    Berlin_13304 <- readRDS("/mnt/Sara/RMS/write/Berlin_13304.rds")
    Rh70 <- readRDS("/mnt/Sara/RMS/write/Rh70.rds")
    Rh73 <- readRDS("/mnt/Sara/RMS/write/Rh73.rds")
    PDX82 <- readRDS("/mnt/Sara/RMS/write/PDX82.rds")
    Rh71 <- readRDS("/mnt/Sara/RMS/write/Rh71.rds")
    Rh74 <- readRDS("/mnt/Sara/RMS/write/Rh74.rds")
    Mast139 <- readRDS("/mnt/Sara/RMS/write/Mast139.rds")
    Berlin_13454 <- readRDS("/mnt/Sara/RMS/write/Berlin_13454.rds")
    Berlin_13933 <- readRDS("/mnt/Sara/RMS/write/Berlin_13933.rds")
    Berlin_13870 <- readRDS("/mnt/Sara/RMS/write/Berlin_13870.rds")
    
    Rh41 <- readRDS("/mnt/Sara/RMS/write/Rh41.rds")

  # Take RNA  
    DefaultAssay(KFR_2) <- "RNA"
    DefaultAssay(Rh4_2) <- "RNA"
    DefaultAssay(RMS_2) <- "RNA"
    DefaultAssay(Rh41) <- "RNA"
    DefaultAssay(PDX29_TS) <- "RNA"
    DefaultAssay(PDX35_TS) <- "RNA"
    DefaultAssay(PDX104_TS) <- "RNA"
    DefaultAssay(Mast118) <- "RNA"
    DefaultAssay(Berlin_13304) <- "RNA"
    DefaultAssay(Rh70) <- "RNA"
    DefaultAssay(Rh73) <- "RNA"
    DefaultAssay(PDX82) <- "RNA"
    DefaultAssay(Rh71) <- "RNA"
    DefaultAssay(Rh74) <- "RNA"
    DefaultAssay(Mast139) <- "RNA"
    DefaultAssay(Berlin_13454) <- "RNA"
    DefaultAssay(Berlin_13870) <- "RNA"
    DefaultAssay(Berlin_13933) <- "RNA"
    
    ## Reduce dimensions of objects
    DietSeurat(KFR_2, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(RMS_2, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Rh4_2, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Rh41, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(PDX104_TS, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(PDX29_TS, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(PDX35_TS, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Mast139, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Mast118, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Rh74, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Rh71, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Rh70, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Rh73, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(PDX82, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Berlin_13933, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Berlin_13454, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Berlin_13304, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    DietSeurat(Berlin_13870, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    
    
    
    # Subset the cells to 1500 cells
    set.seed(111)
    KFR_2_small <- KFR_2[, sample(colnames(KFR_2), size = 1500, replace=F)]
    rm(KFR_2) 
    
    set.seed(111)
    Rh4_2_small <- Rh4_2[, sample(colnames(Rh4_2), size = 1500, replace=F)]
    rm(Rh4_2) 
    
    set.seed(111)
    RMS_2_small <- RMS_2[, sample(colnames(RMS_2), size = 1500, replace=F)]
    rm(RMS_2) 
    
    set.seed(111)
    PDX35_TS_small <- PDX35_TS[, sample(colnames(PDX35_TS), size = 1500, replace=F)]
    rm(PDX35_TS) 
    
    set.seed(111)
    PDX29_TS_small <- PDX29_TS[, sample(colnames(PDX29_TS), size = 1500, replace=F)]
    rm(PDX29_TS) 
    
    set.seed(111)
    PDX104_TS_small <- PDX104_TS[, sample(colnames(PDX104_TS), size = 1500, replace=F)]
    rm(PDX104_TS) 
    
    set.seed(111)
    Rh70_small <- Rh70[, sample(colnames(Rh70), size = 993, replace=F)]
    rm(Rh70) 
    
    set.seed(111)
    Rh73_small <- Rh73[, sample(colnames(Rh73), size = 1500, replace=F)]
    rm(Rh73) 
    
    set.seed(111)
    PDX82_small <- PDX82[, sample(colnames(PDX82), size = 1500, replace=F)]
    rm(PDX82) 
    
    set.seed(111)
    Berlin_13454_small <- Berlin_13454[, sample(colnames(Berlin_13454), size = 1500, replace=F)]
    rm(Berlin_13454) 
    
    set.seed(111)
    Berlin_13933_small <- Berlin_13933[, sample(colnames(Berlin_13933), size = 1500, replace=F)]
    rm(Berlin_13933) 
    
    set.seed(111)
    Berlin_13870_small <- Berlin_13870[, sample(colnames(Berlin_13870), size = 1500, replace=F)]
    rm(Berlin_13870) 
    
    set.seed(111)
    Rh71_small <- Rh71[, sample(colnames(Rh71), size = 1328, replace=F)]
    rm(Rh71) 
    
    set.seed(111)
    Rh74_small <- Rh74[, sample(colnames(Rh74), size = 1500, replace=F)]
    rm(Rh74) 
    
    set.seed(111)
    Mast139_small <- Mast139[, sample(colnames(Mast139), size = 1500, replace=F)]
    rm(Mast139) 
    
    set.seed(111)
    Mast118_small <- Mast118[, sample(colnames(Mast118), size = 1276, replace=F)]
    rm(Mast118) 
    
    set.seed(111)
    Berlin_13304_small <- Berlin_13304[, sample(colnames(Berlin_13304), size = 1865, replace=F)]
    rm(Berlin_13304) 
    
    set.seed(111)
    Rh41_small <- Rh41[, sample(colnames(Rh41), size = 1500, replace=F)]
    rm(Rh41) 
    

  
##################################################################
############ (4) Add metadata ##############
##################################################################
    
# NAME
    PDX104_TS_small@meta.data$name = 'aRMS-1'
    PDX29_TS_small@meta.data$name = 'aRMS-2'
    PDX35_TS_small@meta.data$name = 'aRMS-3'
    Berlin_13304_small@meta.data$name = 'aRMS-5'
    PDX82_small@meta.data$name = 'eRMS-4'
    Berlin_13454_small@meta.data$name = 'eRMS-8.1'
    Berlin_13870_small@meta.data$name = 'eRMS-8.2'
    Berlin_13933_small@meta.data$name = 'eRMS-8.3'
    KFR_2_small@meta.data$name = 'KFR'
    Rh4_2_small@meta.data$name = 'Rh4'
    RMS_2_small@meta.data$name = 'RMS'
    Rh70_small@meta.data$name = 'eRMS-1.1'
    Rh73_small@meta.data$name = 'eRMS-1.2'
    Rh71_small@meta.data$name = 'eRMS-2.1'
    Rh74_small@meta.data$name = 'eRMS-2.2'
    Mast139_small@meta.data$name = 'eRMS-3.2'
    Mast118_small@meta.data$name = 'aRMS-4'
    sm_Langenau_20082@meta.data$name = '20082'
    sm_Langenau_20696@meta.data$name = '20696'
    sm_Langenau_21202@meta.data$name = '21202'
    sm_Langenau_29806@meta.data$name = '29806'
    sm_Langenau_Mast39@meta.data$name = 'Mast39'
    sm_Langenau_MSK72117@meta.data$name = 'MSK72117'
    sm_Langenau_MSK72117_SC@meta.data$name = 'MSK72117_SC'
    sm_Langenau_MSK74711@meta.data$name = 'MSK74711'
    sm_Langenau_MSK82489@meta.data$name = 'MSK82489'
    sm_Langenau_RD@meta.data$name = 'RD'
    sm_Langenau_Mast85_r1@meta.data$name = 'Mast85_r1'
    sm_Langenau_Mast85_r2@meta.data$name = 'Mast85_r2'
    sm_Langenau_Mast85_r2_SC@meta.data$name = 'Mast85_r2_SC'
    sm_Langenau_Mast95@meta.data$name = 'Mast95'
    sm_Langenau_Mast111@meta.data$name = 'Mast111'
    sm_Langenau_Mast139@meta.data$name = 'Mast139'
    sm_Langenau_Mast139_SC@meta.data$name = 'Mast139_SC'
    sm_Langenau_Mast118@meta.data$name = 'Mast118'
    SJRHB000026_R2_snRNAseq_small@meta.data$name = 'SJRHB000026_R2'
    SJRHB000026_X1_scRNAseq_small@meta.data$name = 'SJRHB000026_X1'
    SJRHB000026_R3_snRNAseq_small@meta.data$name = 'SJRHB000026_R3'
    SJRHB000026_X2_scRNAseq_small@meta.data$name = 'SJRHB000026_X2'
    SJRHB010468_D1_snRNAseq_small@meta.data$name = 'SJRHB010468_D1'
    SJRHB010468_X1_scRNAseq_small@meta.data$name = 'SJRHB010468_X1'
    SJRHB010927_D1_snRNAseq_small@meta.data$name = 'SJRHB010927_D1'
    SJRHB010927_X1_scRNAseq_small@meta.data$name = 'SJRHB010927_X1'
    SJRHB010928_R1_snRNAseq_small@meta.data$name = 'SJRHB010928_R1'
    SJRHB010928_X1_scRNAseq_small@meta.data$name = 'SJRHB010928_X1'
    SJRHB011_D_snRNAseq_small@meta.data$name = 'SJRHB011_D'
    SJRHB011_X_scRNAseq_small@meta.data$name = 'SJRHB011_X'
    SJRHB012_R_snRNAseq_small@meta.data$name = 'SJRHB012_R'
    SJRHB012_Y_scRNAseq_small@meta.data$name = 'SJRHB012_Y'
    SJRHB012_S_snRNAseq_small@meta.data$name = 'SJRHB012_S'
    SJRHB012_Z_scRNAseq_small@meta.data$name = 'SJRHB012_Z'
    SJRHB012405_D1_snRNAseq_small@meta.data$name = 'SJRHB012405_D1'
    SJRHB012405_X1_scRNAseq_small@meta.data$name = 'SJRHB012405_X1'
    SJRHB013757_D2_snRNAseq_small@meta.data$name = 'SJRHB013757_D2'
    SJRHB013757_X1_scRNAseq_small@meta.data$name = 'SJRHB013757_X1'
    SJRHB013758_D1_snRNAseq_small@meta.data$name = 'SJRHB013758_D1'
    SJRHB013758_X1_scRNAseq_small@meta.data$name = 'SJRHB013758_X1'
    SJRHB013758_D2_snRNAseq_small@meta.data$name = 'SJRHB013758_D2'
    SJRHB013758_X2_scRNAseq_small@meta.data$name = 'SJRHB013758_X2'
    SJRHB013759_A1_snRNAseq_small@meta.data$name = 'SJRHB013759_A1'
    SJRHB013759_X14_scRNAseq_small@meta.data$name = 'SJRHB013759_X14'
    SJRHB013759_A2_snRNAseq_small@meta.data$name = 'SJRHB013759_A2'
    SJRHB013759_X15_scRNAseq_small@meta.data$name = 'SJRHB013759_X15'
    SJRHB030680_R1_scRNAseq_small@meta.data$name = 'SJRHB030680_R1'
    SJRHB030680_X1_scRNAseq_small@meta.data$name = 'SJRHB030680_X1'
    SJRHB031320_D1_scRNAseq_small@meta.data$name = 'SJRHB031320_D1'
    SJRHB031320_X1_scRNAseq_small@meta.data$name = 'SJRHB031320_X1'
    SJRHB046156_A1_snRNAseq_small@meta.data$name = 'SJRHB046156_A1'
    SJRHB046156_X1_scRNAseq_small@meta.data$name = 'SJRHB046156_X1'
    SJRHB049189_D1_snRNAseq_small@meta.data$name = 'SJRHB049189_D1'
    SJRHB049189_X1_scRNAseq_small@meta.data$name = 'SJRHB049189_X1'
    Rh41_small@meta.data$name = 'Rh41'
  
    
# PATIENT ID
    PDX104_TS_small@meta.data$PatientID = 'aRMS-1'
    PDX29_TS_small@meta.data$PatientID = 'aRMS-2'
    PDX35_TS_small@meta.data$PatientID = 'aRMS-3'
    Berlin_13304_small@meta.data$PatientID = 'aRMS-5'
    PDX82_small@meta.data$PatientID = 'eRMS-4'
    Berlin_13454_small@meta.data$PatientID = 'eRMS-8'
    Berlin_13870_small@meta.data$PatientID = 'eRMS-8'
    Berlin_13933_small@meta.data$PatientID = 'eRMS-8'
    KFR_2_small@meta.data$PatientID = 'KFR'
    Rh4_2_small@meta.data$PatientID = 'Rh4_Rh41'
    RMS_2_small@meta.data$PatientID = 'RMS'
    Rh70_small@meta.data$PatientID = 'SJRHB011'
    Rh73_small@meta.data$PatientID = 'SJRHB011'
    Rh71_small@meta.data$PatientID = 'SJRHB012'
    Rh74_small@meta.data$PatientID = 'SJRHB012'
    Mast139_small@meta.data$PatientID = 'SJRHB013758'
    Mast118_small@meta.data$PatientID = 'SJRHB013759'
    sm_Langenau_20082@meta.data$PatientID = '20082'
    sm_Langenau_20696@meta.data$PatientID = '20696_21202'
    sm_Langenau_21202@meta.data$PatientID = '20696_21202'
    sm_Langenau_29806@meta.data$PatientID = '29806'
    sm_Langenau_Mast39@meta.data$PatientID = 'Mast39'
    sm_Langenau_MSK72117@meta.data$PatientID = 'MSK72117'
    sm_Langenau_MSK72117_SC@meta.data$PatientID = 'MSK72117'
    sm_Langenau_MSK74711@meta.data$PatientID = 'MSK74711'
    sm_Langenau_MSK82489@meta.data$PatientID = 'MSK82489'
    sm_Langenau_RD@meta.data$PatientID = 'RD'
    sm_Langenau_Mast85_r1@meta.data$PatientID = 'MAST85'
    sm_Langenau_Mast85_r2@meta.data$PatientID = 'MAST85'
    sm_Langenau_Mast85_r2_SC@meta.data$PatientID = 'MAST85'
    sm_Langenau_Mast95@meta.data$PatientID = 'SJRHB013757'
    sm_Langenau_Mast111@meta.data$PatientID = 'SJRHB013758'
    sm_Langenau_Mast139@meta.data$PatientID = 'SJRHB013758'
    sm_Langenau_Mast139_SC@meta.data$PatientID = 'SJRHB013758'
    sm_Langenau_Mast118@meta.data$PatientID = 'SJRHB013759'
    SJRHB000026_R2_snRNAseq_small@meta.data$PatientID = 'SJRHB000026'
    SJRHB000026_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB000026'
    SJRHB000026_R3_snRNAseq_small@meta.data$PatientID = 'SJRHB000026'
    SJRHB000026_X2_scRNAseq_small@meta.data$PatientID = 'SJRHB000026'
    SJRHB010468_D1_snRNAseq_small@meta.data$PatientID = 'SJRHB010468'
    SJRHB010468_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB010468'
    SJRHB010927_D1_snRNAseq_small@meta.data$PatientID = 'SJRHB010927'
    SJRHB010927_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB010927'
    SJRHB010928_R1_snRNAseq_small@meta.data$PatientID = 'SJRHB010928'
    SJRHB010928_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB010928'
    SJRHB011_D_snRNAseq_small@meta.data$PatientID = 'SJRHB011'
    SJRHB011_X_scRNAseq_small@meta.data$PatientID = 'SJRHB011'
    SJRHB012_R_snRNAseq_small@meta.data$PatientID = 'SJRHB012'
    SJRHB012_Y_scRNAseq_small@meta.data$PatientID = 'SJRHB012'
    SJRHB012_S_snRNAseq_small@meta.data$PatientID = 'SJRHB012'
    SJRHB012_Z_scRNAseq_small@meta.data$PatientID = 'SJRHB012'
    SJRHB012405_D1_snRNAseq_small@meta.data$PatientID = 'SJRHB012405'
    SJRHB012405_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB012405'
    SJRHB013757_D2_snRNAseq_small@meta.data$PatientID = 'SJRHB013757'
    SJRHB013757_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB013757'
    SJRHB013758_D1_snRNAseq_small@meta.data$PatientID = 'SJRHB013758'
    SJRHB013758_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB013758'
    SJRHB013758_D2_snRNAseq_small@meta.data$PatientID = 'SJRHB013758'
    SJRHB013758_X2_scRNAseq_small@meta.data$PatientID = 'SJRHB013758'
    SJRHB013759_A1_snRNAseq_small@meta.data$PatientID = 'SJRHB013759'
    SJRHB013759_X14_scRNAseq_small@meta.data$PatientID = 'SJRHB013759'
    SJRHB013759_A2_snRNAseq_small@meta.data$PatientID = 'SJRHB013759'
    SJRHB013759_X15_scRNAseq_small@meta.data$PatientID = 'SJRHB013759'
    SJRHB030680_R1_scRNAseq_small@meta.data$PatientID = 'SJRHB030680'
    SJRHB030680_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB030680'
    SJRHB031320_D1_scRNAseq_small@meta.data$PatientID = 'SJRHB031320'
    SJRHB031320_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB031320'
    SJRHB046156_A1_snRNAseq_small@meta.data$PatientID = 'SJRHB046156'
    SJRHB046156_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB046156'
    SJRHB049189_D1_snRNAseq_small@meta.data$PatientID = 'SJRHB049189'
    SJRHB049189_X1_scRNAseq_small@meta.data$PatientID = 'SJRHB049189'
    Rh41_small@meta.data$PatientID = 'Rh4_Rh41'
    
    
# ID    
    PDX104_TS_small@meta.data$id = 'aRMS-1'
    PDX29_TS_small@meta.data$id = 'aRMS-2'
    PDX35_TS_small@meta.data$id = 'aRMS-3'
    Berlin_13304_small@meta.data$id = 'aRMS-5'
    PDX82_small@meta.data$id = 'eRMS-4'
    Berlin_13454_small@meta.data$id = 'eRMS-8.1'
    Berlin_13870_small@meta.data$id = 'eRMS-8.2'
    Berlin_13933_small@meta.data$id = 'eRMS-8.3'
    KFR_2_small@meta.data$id = 'KFR'
    Rh4_2_small@meta.data$id = 'Rh4'
    RMS_2_small@meta.data$id = 'RMS'
    Rh70_small@meta.data$id = 'SJRHB011_D'
    Rh73_small@meta.data$id = 'SJRHB011_E'
    Rh71_small@meta.data$id = 'SJRHB012_R'
    Rh74_small@meta.data$id = 'SJRHB012_S'
    Mast139_small@meta.data$id = 'SJRHB013758_D2'
    Mast118_small@meta.data$id = 'SJRHB013759_A1'
    sm_Langenau_20082@meta.data$id = '20082'
    sm_Langenau_20696@meta.data$id = '20696'
    sm_Langenau_21202@meta.data$id = '21202'
    sm_Langenau_29806@meta.data$id = '29806'
    sm_Langenau_Mast39@meta.data$id = 'Mast39'
    sm_Langenau_MSK72117@meta.data$id = 'MSK72117'
    sm_Langenau_MSK72117_SC@meta.data$id = 'MSK72117'
    sm_Langenau_MSK74711@meta.data$id = 'MSK74711'
    sm_Langenau_MSK82489@meta.data$id = 'MSK82489'
    sm_Langenau_RD@meta.data$id = 'RD'
    sm_Langenau_Mast85_r1@meta.data$id = 'SJRHB000026_R3'
    sm_Langenau_Mast85_r2@meta.data$id = 'SJRHB000026_R3'
    sm_Langenau_Mast85_r2_SC@meta.data$id = 'SJRHB000026_R3'
    sm_Langenau_Mast95@meta.data$id = 'SJRHB013757_D2'
    sm_Langenau_Mast111@meta.data$id = 'SJRHB013758_D1'
    sm_Langenau_Mast139@meta.data$id = 'SJRHB013758_D2'
    sm_Langenau_Mast139_SC@meta.data$id = 'SJRHB013758_D2'
    sm_Langenau_Mast118@meta.data$id = 'SJRHB013759_D1'
    SJRHB000026_R2_snRNAseq_small@meta.data$id = 'SJRHB000026_R2'
    SJRHB000026_X1_scRNAseq_small@meta.data$id = 'SJRHB000026_R2'
    SJRHB000026_R3_snRNAseq_small@meta.data$id = 'SJRHB000026_R3'
    SJRHB000026_X2_scRNAseq_small@meta.data$id = 'SJRHB000026_R3'
    SJRHB010468_D1_snRNAseq_small@meta.data$id = 'SJRHB010468_D1'
    SJRHB010468_X1_scRNAseq_small@meta.data$id = 'SJRHB010468_D1'
    SJRHB010927_D1_snRNAseq_small@meta.data$id = 'SJRHB010927_D1'
    SJRHB010927_X1_scRNAseq_small@meta.data$id = 'SJRHB010927_D1'
    SJRHB010928_R1_snRNAseq_small@meta.data$id = 'SJRHB010928_R1'
    SJRHB010928_X1_scRNAseq_small@meta.data$id = 'SJRHB010928_R1'
    SJRHB011_D_snRNAseq_small@meta.data$id = 'SJRHB011_D'
    SJRHB011_X_scRNAseq_small@meta.data$id = 'SJRHB011_D'
    SJRHB012_R_snRNAseq_small@meta.data$id = 'SJRHB012_R'
    SJRHB012_Y_scRNAseq_small@meta.data$id = 'SJRHB012_R'
    SJRHB012_S_snRNAseq_small@meta.data$id = 'SJRHB012_S'
    SJRHB012_Z_scRNAseq_small@meta.data$id = 'SJRHB012_S'
    SJRHB012405_D1_snRNAseq_small@meta.data$id = 'SJRHB012405_D1'
    SJRHB012405_X1_scRNAseq_small@meta.data$id = 'SJRHB012405_D1'
    SJRHB013757_D2_snRNAseq_small@meta.data$id = 'SJRHB013757_D2'
    SJRHB013757_X1_scRNAseq_small@meta.data$id = 'SJRHB013757_D2'
    SJRHB013758_D1_snRNAseq_small@meta.data$id = 'SJRHB013758_D1'
    SJRHB013758_X1_scRNAseq_small@meta.data$id = 'SJRHB013758_D1'
    SJRHB013758_D2_snRNAseq_small@meta.data$id = 'SJRHB013758_D2'
    SJRHB013758_X2_scRNAseq_small@meta.data$id = 'SJRHB013758_D2'
    SJRHB013759_A1_snRNAseq_small@meta.data$id = 'SJRHB013759_A1'
    SJRHB013759_X14_scRNAseq_small@meta.data$id = 'SJRHB013759_A1'
    SJRHB013759_A2_snRNAseq_small@meta.data$id = 'SJRHB013759_A2'
    SJRHB013759_X15_scRNAseq_small@meta.data$id = 'SJRHB013759_A2'
    SJRHB030680_R1_scRNAseq_small@meta.data$id = 'SJRHB030680_R1'
    SJRHB030680_X1_scRNAseq_small@meta.data$id = 'SJRHB030680_R1'
    SJRHB031320_D1_scRNAseq_small@meta.data$id = 'SJRHB031320_D1'
    SJRHB031320_X1_scRNAseq_small@meta.data$id = 'SJRHB031320_D1'
    SJRHB046156_A1_snRNAseq_small@meta.data$id = 'SJRHB046156_R1'
    SJRHB046156_X1_scRNAseq_small@meta.data$id = 'SJRHB046156_R1'
    SJRHB049189_D1_snRNAseq_small@meta.data$id = 'SJRHB049189_D1'
    SJRHB049189_X1_scRNAseq_small@meta.data$id = 'SJRHB049189_D1'
    Rh41_small@meta.data$id = 'Rh41'
    
  
# ORIGIN      
    PDX104_TS_small@meta.data$origin = 'Danielli et al.'
    PDX29_TS_small@meta.data$origin = 'Danielli et al.'
    PDX35_TS_small@meta.data$origin = 'Danielli et al.'
    Berlin_13304_small@meta.data$origin = 'Danielli et al.'
    PDX82_small@meta.data$origin = 'Danielli et al.'
    Berlin_13454_small@meta.data$origin = 'Danielli et al.'
    Berlin_13870_small@meta.data$origin = 'Danielli et al.'
    Berlin_13933_small@meta.data$origin = 'Danielli et al.'
    KFR_2_small@meta.data$origin = 'Danielli et al.'
    Rh4_2_small@meta.data$origin = 'Danielli et al.'
    RMS_2_small@meta.data$origin = 'Danielli et al.'
    Rh70_small@meta.data$origin = 'Danielli et al.'
    Rh73_small@meta.data$origin = 'Danielli et al.'
    Rh71_small@meta.data$origin = 'Danielli et al.'
    Rh74_small@meta.data$origin = 'Danielli et al.'
    Mast139_small@meta.data$origin = 'Danielli et al.'
    Mast118_small@meta.data$origin = 'Danielli et al.'
    sm_Langenau_20082@meta.data$origin = 'Wei et al.'
    sm_Langenau_20696@meta.data$origin = 'Wei et al.'
    sm_Langenau_21202@meta.data$origin = 'Wei et al.'
    sm_Langenau_29806@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast39@meta.data$origin = 'Wei et al.'
    sm_Langenau_MSK72117@meta.data$origin = 'Wei et al.'
    sm_Langenau_MSK72117_SC@meta.data$origin = 'Wei et al.'
    sm_Langenau_MSK74711@meta.data$origin = 'Wei et al.'
    sm_Langenau_MSK82489@meta.data$origin = 'Wei et al.'
    sm_Langenau_RD@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast85_r1@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast85_r2@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast85_r2_SC@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast95@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast111@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast139@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast139_SC@meta.data$origin = 'Wei et al.'
    sm_Langenau_Mast118@meta.data$origin = 'Wei et al.'
    SJRHB000026_R2_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB000026_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB000026_R3_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB000026_X2_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB010468_D1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB010468_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB010927_D1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB010927_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB010928_R1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB010928_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB011_D_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB011_X_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB012_R_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB012_Y_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB012_S_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB012_Z_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB012405_D1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB012405_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013757_D2_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013757_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013758_D1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013758_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013758_D2_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013758_X2_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013759_A1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013759_X14_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013759_A2_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB013759_X15_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB030680_R1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB030680_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB031320_D1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB031320_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB046156_A1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB046156_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB049189_D1_snRNAseq_small@meta.data$origin = 'Patel et al.'
    SJRHB049189_X1_scRNAseq_small@meta.data$origin = 'Patel et al.'
    Rh41_small@meta.data$origin = 'Weng et al.'
  
    
    
# MODEL
    PDX104_TS_small@meta.data$model = 'Primary culture'
    PDX29_TS_small@meta.data$model = 'Primary culture'
    PDX35_TS_small@meta.data$model = 'Primary culture'
    Berlin_13304_small@meta.data$model = 'Primary culture'
    PDX82_small@meta.data$model = 'Primary culture'
    Berlin_13454_small@meta.data$model = 'Primary culture'
    Berlin_13870_small@meta.data$model = 'Primary culture'
    Berlin_13933_small@meta.data$model = 'Primary culture'
    KFR_2_small@meta.data$model = 'Cell line'
    Rh4_2_small@meta.data$model = 'Cell line'
    RMS_2_small@meta.data$model = 'Cell line'
    Rh70_small@meta.data$model = 'Primary culture'
    Rh73_small@meta.data$model = 'Primary culture'
    Rh71_small@meta.data$model = 'Primary culture'
    Rh74_small@meta.data$model = 'Primary culture'
    Mast139_small@meta.data$model = 'Primary culture'
    Mast118_small@meta.data$model = 'Primary culture'
    sm_Langenau_20082@meta.data$model = 'Patient'
    sm_Langenau_20696@meta.data$model = 'Patient'
    sm_Langenau_21202@meta.data$model = 'Patient'
    sm_Langenau_29806@meta.data$model = 'Patient'
    sm_Langenau_Mast39@meta.data$model = 'O-PDX'
    sm_Langenau_MSK72117@meta.data$model = 'O-PDX'
    sm_Langenau_MSK72117_SC@meta.data$model = 'O-PDX'
    sm_Langenau_MSK74711@meta.data$model = 'O-PDX'
    sm_Langenau_MSK82489@meta.data$model = 'O-PDX'
    sm_Langenau_RD@meta.data$model = 'Cell line'
    sm_Langenau_Mast85_r1@meta.data$model = 'O-PDX'
    sm_Langenau_Mast85_r2@meta.data$model = 'O-PDX'
    sm_Langenau_Mast85_r2_SC@meta.data$model = 'O-PDX'
    sm_Langenau_Mast95@meta.data$model = 'O-PDX'
    sm_Langenau_Mast111@meta.data$model = 'O-PDX'
    sm_Langenau_Mast139@meta.data$model = 'O-PDX'
    sm_Langenau_Mast139_SC@meta.data$model = 'O-PDX'
    sm_Langenau_Mast118@meta.data$model = 'O-PDX'
    SJRHB000026_R2_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB000026_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB000026_R3_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB000026_X2_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB010468_D1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB010468_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB010927_D1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB010927_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB010928_R1_snRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB010928_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB011_D_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB011_X_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB012_R_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB012_Y_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB012_S_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB012_Z_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB012405_D1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB012405_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB013757_D2_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB013757_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB013758_D1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB013758_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB013758_D2_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB013758_X2_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB013759_A1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB013759_X14_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB013759_A2_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB013759_X15_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB030680_R1_scRNAseq_small@meta.data$model = 'Patient'
    SJRHB030680_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB031320_D1_scRNAseq_small@meta.data$model = 'Patient'
    SJRHB031320_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB046156_A1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB046156_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    SJRHB049189_D1_snRNAseq_small@meta.data$model = 'Patient'
    SJRHB049189_X1_scRNAseq_small@meta.data$model = 'O-PDX'
    Rh41_small@meta.data$model = 'Cell line'
    
    
# SEQUENCING  
    PDX104_TS_small@meta.data$sequencing = 'scRNAseq'
    PDX29_TS_small@meta.data$sequencing = 'scRNAseq'
    PDX35_TS_small@meta.data$sequencing = 'scRNAseq'
    Berlin_13304_small@meta.data$sequencing = 'scRNAseq'
    PDX82_small@meta.data$sequencing = 'scRNAseq'
    Berlin_13454_small@meta.data$sequencing = 'scRNAseq'
    Berlin_13870_small@meta.data$sequencing = 'scRNAseq'
    Berlin_13933_small@meta.data$sequencing = 'scRNAseq'
    KFR_2_small@meta.data$sequencing = 'scRNAseq'
    Rh4_2_small@meta.data$sequencing = 'scRNAseq'
    RMS_2_small@meta.data$sequencing = 'scRNAseq'
    Rh70_small@meta.data$sequencing = 'scRNAseq'
    Rh73_small@meta.data$sequencing = 'scRNAseq'
    Rh71_small@meta.data$sequencing = 'scRNAseq'
    Rh74_small@meta.data$sequencing = 'scRNAseq'
    Mast139_small@meta.data$sequencing = 'scRNAseq'
    Mast118_small@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_20082@meta.data$sequencing = 'snRNAseq'
    sm_Langenau_20696@meta.data$sequencing = 'snRNAseq'
    sm_Langenau_21202@meta.data$sequencing = 'snRNAseq'
    sm_Langenau_29806@meta.data$sequencing = 'snRNAseq'
    sm_Langenau_Mast39@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_MSK72117@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_MSK72117_SC@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_MSK74711@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_MSK82489@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_RD@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast85_r1@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast85_r2@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast85_r2_SC@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast95@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast111@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast139@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast139_SC@meta.data$sequencing = 'scRNAseq'
    sm_Langenau_Mast118@meta.data$sequencing = 'scRNAseq'
    SJRHB000026_R2_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB000026_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB000026_R3_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB000026_X2_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB010468_D1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB010468_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB010927_D1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB010927_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB010928_R1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB010928_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB011_D_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB011_X_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB012_R_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB012_Y_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB012_S_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB012_Z_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB012405_D1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB012405_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB013757_D2_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB013757_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB013758_D1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB013758_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB013758_D2_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB013758_X2_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB013759_A1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB013759_X14_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB013759_A2_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB013759_X15_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB030680_R1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB030680_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB031320_D1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB031320_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB046156_A1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB046156_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    SJRHB049189_D1_snRNAseq_small@meta.data$sequencing = 'snRNAseq'
    SJRHB049189_X1_scRNAseq_small@meta.data$sequencing = 'scRNAseq'
    Rh41_small@meta.data$sequencing = 'scRNAseq'
    
  
# SEX   
    PDX104_TS_small@meta.data$sex = 'F'
    PDX29_TS_small@meta.data$sex = 'F'
    PDX35_TS_small@meta.data$sex = 'M'
    Berlin_13304_small@meta.data$sex = 'NA'
    PDX82_small@meta.data$sex = 'F'
    Berlin_13454_small@meta.data$sex = 'F'
    Berlin_13870_small@meta.data$sex = 'F'
    Berlin_13933_small@meta.data$sex = 'F'
    KFR_2_small@meta.data$sex = 'F'
    Rh4_2_small@meta.data$sex = 'F'
    RMS_2_small@meta.data$sex = 'NA'
    Rh70_small@meta.data$sex = 'M'
    Rh73_small@meta.data$sex = 'M'
    Rh71_small@meta.data$sex = 'M'
    Rh74_small@meta.data$sex = 'M'
    Mast139_small@meta.data$sex = 'F'
    Mast118_small@meta.data$sex = 'M'
    sm_Langenau_20082@meta.data$sex = 'F'
    sm_Langenau_20696@meta.data$sex = 'M'
    sm_Langenau_21202@meta.data$sex = 'M'
    sm_Langenau_29806@meta.data$sex = 'M'
    sm_Langenau_Mast39@meta.data$sex = 'F'
    sm_Langenau_MSK72117@meta.data$sex = 'F'
    sm_Langenau_MSK72117_SC@meta.data$sex = 'F'
    sm_Langenau_MSK74711@meta.data$sex = 'F'
    sm_Langenau_MSK82489@meta.data$sex = 'F'
    sm_Langenau_RD@meta.data$sex = 'NA'
    sm_Langenau_Mast85_r1@meta.data$sex = 'F'
    sm_Langenau_Mast85_r2@meta.data$sex = 'F'
    sm_Langenau_Mast85_r2_SC@meta.data$sex = 'F'
    sm_Langenau_Mast95@meta.data$sex = 'F'
    sm_Langenau_Mast111@meta.data$sex = 'F'
    sm_Langenau_Mast139@meta.data$sex = 'F'
    sm_Langenau_Mast139_SC@meta.data$sex = 'F'
    sm_Langenau_Mast118@meta.data$sex = 'M'
    SJRHB000026_R2_snRNAseq_small@meta.data$sex = 'F'
    SJRHB000026_X1_scRNAseq_small@meta.data$sex = 'F'
    SJRHB000026_R3_snRNAseq_small@meta.data$sex = 'F'
    SJRHB000026_X2_scRNAseq_small@meta.data$sex = 'F'
    SJRHB010468_D1_snRNAseq_small@meta.data$sex = 'M'
    SJRHB010468_X1_scRNAseq_small@meta.data$sex = 'M'
    SJRHB010927_D1_snRNAseq_small@meta.data$sex = 'F'
    SJRHB010927_X1_scRNAseq_small@meta.data$sex = 'F'
    SJRHB010928_R1_snRNAseq_small@meta.data$sex = 'M'
    SJRHB010928_X1_scRNAseq_small@meta.data$sex = 'M'
    SJRHB011_D_snRNAseq_small@meta.data$sex = 'M'
    SJRHB011_X_scRNAseq_small@meta.data$sex = 'M'
    SJRHB012_R_snRNAseq_small@meta.data$sex = 'M'
    SJRHB012_Y_scRNAseq_small@meta.data$sex = 'M'
    SJRHB012_S_snRNAseq_small@meta.data$sex = 'M'
    SJRHB012_Z_scRNAseq_small@meta.data$sex = 'M'
    SJRHB012405_D1_snRNAseq_small@meta.data$sex = 'F'
    SJRHB012405_X1_scRNAseq_small@meta.data$sex = 'F'
    SJRHB013757_D2_snRNAseq_small@meta.data$sex = 'F'
    SJRHB013757_X1_scRNAseq_small@meta.data$sex = 'F'
    SJRHB013758_D1_snRNAseq_small@meta.data$sex = 'F'
    SJRHB013758_X1_scRNAseq_small@meta.data$sex = 'F'
    SJRHB013758_D2_snRNAseq_small@meta.data$sex = 'F'
    SJRHB013758_X2_scRNAseq_small@meta.data$sex = 'F'
    SJRHB013759_A1_snRNAseq_small@meta.data$sex = 'M'
    SJRHB013759_X14_scRNAseq_small@meta.data$sex = 'M'
    SJRHB013759_A2_snRNAseq_small@meta.data$sex = 'M'
    SJRHB013759_X15_scRNAseq_small@meta.data$sex = 'M'
    SJRHB030680_R1_scRNAseq_small@meta.data$sex = 'M'
    SJRHB030680_X1_scRNAseq_small@meta.data$sex = 'M'
    SJRHB031320_D1_scRNAseq_small@meta.data$sex = 'M'
    SJRHB031320_X1_scRNAseq_small@meta.data$sex = 'M'
    SJRHB046156_A1_snRNAseq_small@meta.data$sex = 'F'
    SJRHB046156_X1_scRNAseq_small@meta.data$sex = 'F'
    SJRHB049189_D1_snRNAseq_small@meta.data$sex = 'M'
    SJRHB049189_X1_scRNAseq_small@meta.data$sex = 'M'
    Rh41_small@meta.data$sex = 'F'
 
    
    
# AGE 
    PDX104_TS_small@meta.data$age = '7'
    PDX29_TS_small@meta.data$age = '14'
    PDX35_TS_small@meta.data$age = '13'
    Berlin_13304_small@meta.data$age = 'NA'
    PDX82_small@meta.data$age = '14'
    Berlin_13454_small@meta.data$age = '<1'
    Berlin_13870_small@meta.data$age = '<1'
    Berlin_13933_small@meta.data$age = '1'
    KFR_2_small@meta.data$age = '13'
    Rh4_2_small@meta.data$age = '7'
    RMS_2_small@meta.data$age = '4'
    Rh70_small@meta.data$age = '4'
    Rh73_small@meta.data$age = '5'
    Rh71_small@meta.data$age = '17'
    Rh74_small@meta.data$age = '18'
    Mast139_small@meta.data$age = '5'
    Mast118_small@meta.data$age = '19'
    sm_Langenau_20082@meta.data$age = '17'
    sm_Langenau_20696@meta.data$age = '<1'
    sm_Langenau_21202@meta.data$age = '1'
    sm_Langenau_29806@meta.data$age = '15'
    sm_Langenau_Mast39@meta.data$age = '4'
    sm_Langenau_MSK72117@meta.data$age = '1'
    sm_Langenau_MSK72117_SC@meta.data$age = '1'
    sm_Langenau_MSK74711@meta.data$age = '44'
    sm_Langenau_MSK82489@meta.data$age = '12'
    sm_Langenau_RD@meta.data$age = 'NA'
    sm_Langenau_Mast85_r1@meta.data$age = '5'
    sm_Langenau_Mast85_r2@meta.data$age = '5'
    sm_Langenau_Mast85_r2_SC@meta.data$age = '5'
    sm_Langenau_Mast95@meta.data$age = '3'
    sm_Langenau_Mast111@meta.data$age = '4'
    sm_Langenau_Mast139@meta.data$age = '4'
    sm_Langenau_Mast139_SC@meta.data$age = '4'
    sm_Langenau_Mast118@meta.data$age = '19'
    SJRHB000026_R2_snRNAseq_small@meta.data$age = '4'
    SJRHB000026_X1_scRNAseq_small@meta.data$age = '4'
    SJRHB000026_R3_snRNAseq_small@meta.data$age = '5'
    SJRHB000026_X2_scRNAseq_small@meta.data$age = '5'
    SJRHB010468_D1_snRNAseq_small@meta.data$age = '1'
    SJRHB010468_X1_scRNAseq_small@meta.data$age = '1'
    SJRHB010927_D1_snRNAseq_small@meta.data$age = '5'
    SJRHB010927_X1_scRNAseq_small@meta.data$age = '5'
    SJRHB010928_R1_snRNAseq_small@meta.data$age = '9'
    SJRHB010928_X1_scRNAseq_small@meta.data$age = '9'
    SJRHB011_D_snRNAseq_small@meta.data$age = '5'
    SJRHB011_X_scRNAseq_small@meta.data$age = '5'
    SJRHB012_R_snRNAseq_small@meta.data$age = '18'
    SJRHB012_Y_scRNAseq_small@meta.data$age = '18'
    SJRHB012_S_snRNAseq_small@meta.data$age = '18'
    SJRHB012_Z_scRNAseq_small@meta.data$age = '18'
    SJRHB012405_D1_snRNAseq_small@meta.data$age = '8'
    SJRHB012405_X1_scRNAseq_small@meta.data$age = '8'
    SJRHB013757_D2_snRNAseq_small@meta.data$age = '3'
    SJRHB013757_X1_scRNAseq_small@meta.data$age = '3'
    SJRHB013758_D1_snRNAseq_small@meta.data$age = '4'
    SJRHB013758_X1_scRNAseq_small@meta.data$age = '4'
    SJRHB013758_D2_snRNAseq_small@meta.data$age = '5'
    SJRHB013758_X2_scRNAseq_small@meta.data$age = '5'
    SJRHB013759_A1_snRNAseq_small@meta.data$age = '19'
    SJRHB013759_X14_scRNAseq_small@meta.data$age = '19'
    SJRHB013759_A2_snRNAseq_small@meta.data$age = '19'
    SJRHB013759_X15_scRNAseq_small@meta.data$age = '19'
    SJRHB030680_R1_scRNAseq_small@meta.data$age = '1'
    SJRHB030680_X1_scRNAseq_small@meta.data$age = '1'
    SJRHB031320_D1_scRNAseq_small@meta.data$age = '17'
    SJRHB031320_X1_scRNAseq_small@meta.data$age = '17'
    SJRHB046156_A1_snRNAseq_small@meta.data$age = '16'
    SJRHB046156_X1_scRNAseq_small@meta.data$age = '16'
    SJRHB049189_D1_snRNAseq_small@meta.data$age = '<1'
    SJRHB049189_X1_scRNAseq_small@meta.data$age = '<1'
    Rh41_small@meta.data$age = '7'
      
  
    
# SUBTYPE   
    PDX104_TS_small@meta.data$subtype = 'aRMS'
    PDX29_TS_small@meta.data$subtype = 'aRMS'
    PDX35_TS_small@meta.data$subtype = 'aRMS'
    Berlin_13304_small@meta.data$subtype = 'aRMS'
    PDX82_small@meta.data$subtype = 'eRMS'
    Berlin_13454_small@meta.data$subtype = 'eRMS'
    Berlin_13870_small@meta.data$subtype = 'eRMS'
    Berlin_13933_small@meta.data$subtype = 'eRMS'
    KFR_2_small@meta.data$subtype = 'aRMS'
    Rh4_2_small@meta.data$subtype = 'aRMS'
    RMS_2_small@meta.data$subtype = 'aRMS'
    Rh70_small@meta.data$subtype = 'eRMS'
    Rh73_small@meta.data$subtype = 'eRMS'
    Rh71_small@meta.data$subtype = 'eRMS'
    Rh74_small@meta.data$subtype = 'eRMS'
    Mast139_small@meta.data$subtype = 'eRMS'
    Mast118_small@meta.data$subtype = 'aRMS'
    sm_Langenau_20082@meta.data$subtype = 'aRMS'
    sm_Langenau_20696@meta.data$subtype = 'eRMS'
    sm_Langenau_21202@meta.data$subtype = 'eRMS'
    sm_Langenau_29806@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast39@meta.data$subtype = 'eRMS'
    sm_Langenau_MSK72117@meta.data$subtype = 'aRMS'
    sm_Langenau_MSK72117_SC@meta.data$subtype = 'aRMS'
    sm_Langenau_MSK74711@meta.data$subtype = 'eRMS'
    sm_Langenau_MSK82489@meta.data$subtype = 'aRMS'
    sm_Langenau_RD@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast85_r1@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast85_r2@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast85_r2_SC@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast95@meta.data$subtype = 'aRMS'
    sm_Langenau_Mast111@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast139@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast139_SC@meta.data$subtype = 'eRMS'
    sm_Langenau_Mast118@meta.data$subtype = 'aRMS'
    SJRHB000026_R2_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB000026_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB000026_R3_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB000026_X2_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB010468_D1_snRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB010468_X1_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB010927_D1_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB010927_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB010928_R1_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB010928_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB011_D_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB011_X_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB012_R_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB012_Y_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB012_S_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB012_Z_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB012405_D1_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB012405_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB013757_D2_snRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB013757_X1_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB013758_D1_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB013758_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB013758_D2_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB013758_X2_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB013759_A1_snRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB013759_X14_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB013759_A2_snRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB013759_X15_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB030680_R1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB030680_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB031320_D1_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB031320_X1_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB046156_A1_snRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB046156_X1_scRNAseq_small@meta.data$subtype = 'aRMS'
    SJRHB049189_D1_snRNAseq_small@meta.data$subtype = 'eRMS'
    SJRHB049189_X1_scRNAseq_small@meta.data$subtype = 'eRMS'
    Rh41_small@meta.data$subtype = 'aRMS'
    
    
    
# FUSION     
    PDX104_TS_small@meta.data$fusion= 'PAX3::FOXO1'
    PDX29_TS_small@meta.data$fusion= 'PAX3::FOXO1'
    PDX35_TS_small@meta.data$fusion= 'PAX3::FOXO1'
    Berlin_13304_small@meta.data$fusion= 'PAX3::FOXO1'
    PDX82_small@meta.data$fusion= 'FN-RMS'
    Berlin_13454_small@meta.data$fusion= 'FN-RMS'
    Berlin_13870_small@meta.data$fusion= 'FN-RMS'
    Berlin_13933_small@meta.data$fusion= 'FN-RMS'
    KFR_2_small@meta.data$fusion= 'PAX3::FOXO1'
    Rh4_2_small@meta.data$fusion= 'PAX3::FOXO1'
    RMS_2_small@meta.data$fusion= 'PAX3::FOXO1'
    Rh70_small@meta.data$fusion= 'FN-RMS'
    Rh73_small@meta.data$fusion= 'FN-RMS'
    Rh71_small@meta.data$fusion= 'FN-RMS'
    Rh74_small@meta.data$fusion= 'FN-RMS'
    Mast139_small@meta.data$fusion= 'FN-RMS'
    Mast118_small@meta.data$fusion= 'PAX3::FOXO1'
    sm_Langenau_20082@meta.data$fusion= 'FN-RMS'
    sm_Langenau_20696@meta.data$fusion= 'FN-RMS'
    sm_Langenau_21202@meta.data$fusion= 'FN-RMS'
    sm_Langenau_29806@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast39@meta.data$fusion= 'FN-RMS'
    sm_Langenau_MSK72117@meta.data$fusion= 'PAX7::FOXO1'
    sm_Langenau_MSK72117_SC@meta.data$fusion= 'PAX7::FOXO1'
    sm_Langenau_MSK74711@meta.data$fusion= 'MYOD1'
    sm_Langenau_MSK82489@meta.data$fusion= 'PAX3::FOXO1'
    sm_Langenau_RD@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast85_r1@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast85_r2@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast85_r2_SC@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast95@meta.data$fusion= 'PAX7::FOXO1'
    sm_Langenau_Mast111@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast139@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast139_SC@meta.data$fusion= 'FN-RMS'
    sm_Langenau_Mast118@meta.data$fusion= 'PAX3::FOXO1'
    SJRHB000026_R2_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB000026_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB000026_R3_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB000026_X2_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB010468_D1_snRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB010468_X1_scRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB010927_D1_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB010927_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB010928_R1_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB010928_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB011_D_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB011_X_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB012_R_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB012_Y_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB012_S_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB012_Z_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB012405_D1_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB012405_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB013757_D2_snRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB013757_X1_scRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB013758_D1_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB013758_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB013758_D2_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB013758_X2_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB013759_A1_snRNAseq_small@meta.data$fusion= 'PAX3::FOXO1'
    SJRHB013759_X14_scRNAseq_small@meta.data$fusion= 'PAX3::FOXO1'
    SJRHB013759_A2_snRNAseq_small@meta.data$fusion= 'PAX3::FOXO1'
    SJRHB013759_X15_scRNAseq_small@meta.data$fusion= 'PAX3::FOXO1'
    SJRHB030680_R1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB030680_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB031320_D1_scRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB031320_X1_scRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB046156_A1_snRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB046156_X1_scRNAseq_small@meta.data$fusion= 'PAX7::FOXO1'
    SJRHB049189_D1_snRNAseq_small@meta.data$fusion= 'FN-RMS'
    SJRHB049189_X1_scRNAseq_small@meta.data$fusion= 'FN-RMS'
    Rh41_small@meta.data$fusion= 'PAX3::FOXO1'
    
# SITE  
    PDX104_TS_small@meta.data$site = 'Primary'
    PDX29_TS_small@meta.data$site = 'Primary'
    PDX35_TS_small@meta.data$site = 'Metastasis'
    Berlin_13304_small@meta.data$site = 'NA'
    PDX82_small@meta.data$site = 'Primary'
    Berlin_13454_small@meta.data$site = 'Primary'
    Berlin_13870_small@meta.data$site = 'Primary'
    Berlin_13933_small@meta.data$site = 'Primary'
    KFR_2_small@meta.data$site = 'NA'
    Rh4_2_small@meta.data$site = 'Metastasis'
    RMS_2_small@meta.data$site = 'NA'
    Rh70_small@meta.data$site = 'Primary'
    Rh73_small@meta.data$site = 'Primary'
    Rh71_small@meta.data$site = 'Primary'
    Rh74_small@meta.data$site = 'Primary'
    Mast139_small@meta.data$site = 'Primary'
    Mast118_small@meta.data$site = 'Metastasis'
    sm_Langenau_20082@meta.data$site = 'Metastasis'
    sm_Langenau_20696@meta.data$site = 'Primary'
    sm_Langenau_21202@meta.data$site = 'NA'
    sm_Langenau_29806@meta.data$site = 'NA'
    sm_Langenau_Mast39@meta.data$site = 'Primary'
    sm_Langenau_MSK72117@meta.data$site = 'NA'
    sm_Langenau_MSK72117_SC@meta.data$site = 'NA'
    sm_Langenau_MSK74711@meta.data$site = 'NA'
    sm_Langenau_MSK82489@meta.data$site = 'Metastasis'
    sm_Langenau_RD@meta.data$site = 'NA'
    sm_Langenau_Mast85_r1@meta.data$site = 'Metastasis'
    sm_Langenau_Mast85_r2@meta.data$site = 'Metastasis'
    sm_Langenau_Mast85_r2_SC@meta.data$site = 'Metastasis'
    sm_Langenau_Mast95@meta.data$site = 'Primary'
    sm_Langenau_Mast111@meta.data$site = 'Primary'
    sm_Langenau_Mast139@meta.data$site = 'Primary'
    sm_Langenau_Mast139_SC@meta.data$site = 'Primary'
    sm_Langenau_Mast118@meta.data$site = 'Metastasis'
    SJRHB000026_R2_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB000026_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB000026_R3_snRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB000026_X2_scRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB010468_D1_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB010468_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB010927_D1_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB010927_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB010928_R1_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB010928_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB011_D_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB011_X_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB012_R_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB012_Y_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB012_S_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB012_Z_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB012405_D1_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB012405_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB013757_D2_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB013757_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB013758_D1_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB013758_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB013758_D2_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB013758_X2_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB013759_A1_snRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB013759_X14_scRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB013759_A2_snRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB013759_X15_scRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB030680_R1_scRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB030680_X1_scRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB031320_D1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB031320_X1_scRNAseq_small@meta.data$site = 'Primary'
    SJRHB046156_A1_snRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB046156_X1_scRNAseq_small@meta.data$site = 'Metastasis'
    SJRHB049189_D1_snRNAseq_small@meta.data$site = 'Primary'
    SJRHB049189_X1_scRNAseq_small@meta.data$site = 'Primary'
    Rh41_small@meta.data$site = 'Metastasis'
    
   
# STATUS
    PDX104_TS_small@meta.data$status = 'Recurrent'
    PDX29_TS_small@meta.data$status = 'Recurrent'
    PDX35_TS_small@meta.data$status = 'Recurrent'
    Berlin_13304_small@meta.data$status = 'NA'
    PDX82_small@meta.data$status = 'Diagnosis'
    Berlin_13454_small@meta.data$status = 'Diagnosis'
    Berlin_13870_small@meta.data$status = 'Recurrent'
    Berlin_13933_small@meta.data$status = 'Recurrent'
    KFR_2_small@meta.data$status = 'NA'
    Rh4_2_small@meta.data$status = 'Recurrent'
    RMS_2_small@meta.data$status = 'NA'
    Rh70_small@meta.data$status = 'Diagnosis'
    Rh73_small@meta.data$status = 'Recurrent'
    Rh71_small@meta.data$status = 'Diagnosis'
    Rh74_small@meta.data$status = 'Recurrent'
    Mast139_small@meta.data$status = 'Recurrent'
    Mast118_small@meta.data$status = 'Recurrent'
    sm_Langenau_20082@meta.data$status = 'NA'
    sm_Langenau_20696@meta.data$status = 'Diagnosis'
    sm_Langenau_21202@meta.data$status = 'Recurrent'
    sm_Langenau_29806@meta.data$status = 'Recurrent'
    sm_Langenau_Mast39@meta.data$status = 'Recurrent'
    sm_Langenau_MSK72117@meta.data$status = 'NA'
    sm_Langenau_MSK72117_SC@meta.data$status = 'NA'
    sm_Langenau_MSK74711@meta.data$status = 'NA'
    sm_Langenau_MSK82489@meta.data$status = 'Recurrent'
    sm_Langenau_RD@meta.data$status = 'NA'
    sm_Langenau_Mast85_r1@meta.data$status = 'Recurrent'
    sm_Langenau_Mast85_r2@meta.data$status = 'Recurrent'
    sm_Langenau_Mast85_r2_SC@meta.data$status = 'Recurrent'
    sm_Langenau_Mast95@meta.data$status = 'Diagnosis'
    sm_Langenau_Mast111@meta.data$status = 'Diagnosis'
    sm_Langenau_Mast139@meta.data$status = 'Recurrent'
    sm_Langenau_Mast139_SC@meta.data$status = 'Recurrent'
    sm_Langenau_Mast118@meta.data$status = 'Recurrent'
    SJRHB000026_R2_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB000026_X1_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB000026_R3_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB000026_X2_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB010468_D1_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB010468_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB010927_D1_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB010927_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB010928_R1_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB010928_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB011_D_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB011_X_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB012_R_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB012_Y_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB012_S_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB012_Z_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB012405_D1_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB012405_X1_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB013757_D2_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB013757_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB013758_D1_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB013758_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB013758_D2_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB013758_X2_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB013759_A1_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB013759_X14_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB013759_A2_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB013759_X15_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB030680_R1_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB030680_X1_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB031320_D1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB031320_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB046156_A1_snRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB046156_X1_scRNAseq_small@meta.data$status = 'Recurrent'
    SJRHB049189_D1_snRNAseq_small@meta.data$status = 'Diagnosis'
    SJRHB049189_X1_scRNAseq_small@meta.data$status = 'Diagnosis'
    Rh41_small@meta.data$status = 'Recurrent'
    
    
# LOCATION
    PDX104_TS_small@meta.data$location = 'Tibia'
    PDX29_TS_small@meta.data$location = 'Paravertebral'
    PDX35_TS_small@meta.data$location = 'mediastinum'
    Berlin_13304_small@meta.data$location = 'NA'
    PDX82_small@meta.data$location = 'Oral cavity'
    Berlin_13454_small@meta.data$location = 'Head'
    Berlin_13870_small@meta.data$location = 'Head'
    Berlin_13933_small@meta.data$location = 'Head'
    KFR_2_small@meta.data$location = 'NA'
    Rh4_2_small@meta.data$location = 'Lung'
    RMS_2_small@meta.data$location = 'NA'
    Rh70_small@meta.data$location = 'Skull'
    Rh73_small@meta.data$location = 'Neck'
    Rh71_small@meta.data$location = 'Prostate'
    Rh74_small@meta.data$location = 'Prostate'
    Mast139_small@meta.data$location = 'Abdomen'
    Mast118_small@meta.data$location = 'Inguinal'
    sm_Langenau_20082@meta.data$location = 'NA'
    sm_Langenau_20696@meta.data$location = 'Lung'
    sm_Langenau_21202@meta.data$location = 'Lung'
    sm_Langenau_29806@meta.data$location = 'NA'
    sm_Langenau_Mast39@meta.data$location = 'Pelvis'
    sm_Langenau_MSK72117@meta.data$location = 'Paraspinal'
    sm_Langenau_MSK72117_SC@meta.data$location = 'Paraspinal'
    sm_Langenau_MSK74711@meta.data$location = 'NA'
    sm_Langenau_MSK82489@meta.data$location = 'Inguinal'
    sm_Langenau_RD@meta.data$location = 'NA'
    sm_Langenau_Mast85_r1@meta.data$location = 'Stomach'
    sm_Langenau_Mast85_r2@meta.data$location = 'Stomach'
    sm_Langenau_Mast85_r2_SC@meta.data$location = 'Stomach'
    sm_Langenau_Mast95@meta.data$location = 'Calf'
    sm_Langenau_Mast111@meta.data$location = 'Abdomen-pelvis'
    sm_Langenau_Mast139@meta.data$location = 'Abdomen-pelvis'
    sm_Langenau_Mast139_SC@meta.data$location = 'Abdomen-pelvis'
    sm_Langenau_Mast118@meta.data$location = 'Inguinal'
    SJRHB000026_R2_snRNAseq_small@meta.data$location = 'Pelvis'
    SJRHB000026_X1_scRNAseq_small@meta.data$location = 'Pelvis'
    SJRHB000026_R3_snRNAseq_small@meta.data$location = 'Pelvis'
    SJRHB000026_X2_scRNAseq_small@meta.data$location = 'Pelvis'
    SJRHB010468_D1_snRNAseq_small@meta.data$location = 'Thigh'
    SJRHB010468_X1_scRNAseq_small@meta.data$location = 'Thigh'
    SJRHB010927_D1_snRNAseq_small@meta.data$location = 'Parapharyngeal mass'
    SJRHB010927_X1_scRNAseq_small@meta.data$location = 'Parapharyngeal mass'
    SJRHB010928_R1_snRNAseq_small@meta.data$location = 'Prostate-bladder'
    SJRHB010928_X1_scRNAseq_small@meta.data$location = 'Prostate-bladder'
    SJRHB011_D_snRNAseq_small@meta.data$location = 'Neck '
    SJRHB011_X_scRNAseq_small@meta.data$location = 'Neck '
    SJRHB012_R_snRNAseq_small@meta.data$location = 'Prostate'
    SJRHB012_Y_scRNAseq_small@meta.data$location = 'Prostate'
    SJRHB012_S_snRNAseq_small@meta.data$location = 'Prostate'
    SJRHB012_Z_scRNAseq_small@meta.data$location = 'Prostate'
    SJRHB012405_D1_snRNAseq_small@meta.data$location = 'Abdomen-pelvis'
    SJRHB012405_X1_scRNAseq_small@meta.data$location = 'Abdomen-pelvis'
    SJRHB013757_D2_snRNAseq_small@meta.data$location = 'Calf'
    SJRHB013757_X1_scRNAseq_small@meta.data$location = 'Calf'
    SJRHB013758_D1_snRNAseq_small@meta.data$location = 'Abdomen-pelvis'
    SJRHB013758_X1_scRNAseq_small@meta.data$location = 'Abdomen-pelvis'
    SJRHB013758_D2_snRNAseq_small@meta.data$location = 'Abdomen-pelvis'
    SJRHB013758_X2_scRNAseq_small@meta.data$location = 'Abdomen-pelvis'
    SJRHB013759_A1_snRNAseq_small@meta.data$location = 'Pelvis '
    SJRHB013759_X14_scRNAseq_small@meta.data$location = 'Pelvis '
    SJRHB013759_A2_snRNAseq_small@meta.data$location = 'Pelvis '
    SJRHB013759_X15_scRNAseq_small@meta.data$location = 'Pelvis '
    SJRHB030680_R1_scRNAseq_small@meta.data$location = 'Prostate'
    SJRHB030680_X1_scRNAseq_small@meta.data$location = 'Prostate'
    SJRHB031320_D1_scRNAseq_small@meta.data$location = 'Abdominal'
    SJRHB031320_X1_scRNAseq_small@meta.data$location = 'Abdominal'
    SJRHB046156_A1_snRNAseq_small@meta.data$location = 'Forearm'
    SJRHB046156_X1_scRNAseq_small@meta.data$location = 'Forearm'
    SJRHB049189_D1_snRNAseq_small@meta.data$location = 'Prostate'
    SJRHB049189_X1_scRNAseq_small@meta.data$location = 'Prostate'
    Rh41_small@meta.data$location = 'Lung'
  
    
    
# TREATMENT
    PDX104_TS_small@meta.data$treatment = 'Yes'
    PDX29_TS_small@meta.data$treatment = 'Yes'
    PDX35_TS_small@meta.data$treatment = 'Yes'
    Berlin_13304_small@meta.data$treatment = 'NA'
    PDX82_small@meta.data$treatment = 'NA'
    Berlin_13454_small@meta.data$treatment = 'No'
    Berlin_13870_small@meta.data$treatment = 'Yes'
    Berlin_13933_small@meta.data$treatment = 'Yes'
    KFR_2_small@meta.data$treatment = 'NA'
    Rh4_2_small@meta.data$treatment = 'NA'
    RMS_2_small@meta.data$treatment = 'NA'
    Rh70_small@meta.data$treatment = 'Yes'
    Rh73_small@meta.data$treatment = 'No'
    Rh71_small@meta.data$treatment = 'No'
    Rh74_small@meta.data$treatment = 'Yes'
    Mast139_small@meta.data$treatment = 'Yes'
    Mast118_small@meta.data$treatment = 'No'
    sm_Langenau_20082@meta.data$treatment = 'NA'
    sm_Langenau_20696@meta.data$treatment = 'NA'
    sm_Langenau_21202@meta.data$treatment = 'NA'
    sm_Langenau_29806@meta.data$treatment = 'NA'
    sm_Langenau_Mast39@meta.data$treatment = 'NA'
    sm_Langenau_MSK72117@meta.data$treatment = 'NA'
    sm_Langenau_MSK72117_SC@meta.data$treatment = 'NA'
    sm_Langenau_MSK74711@meta.data$treatment = 'NA'
    sm_Langenau_MSK82489@meta.data$treatment = 'NA'
    sm_Langenau_RD@meta.data$treatment = 'NA'
    sm_Langenau_Mast85_r1@meta.data$treatment = 'NA'
    sm_Langenau_Mast85_r2@meta.data$treatment = 'NA'
    sm_Langenau_Mast85_r2_SC@meta.data$treatment = 'NA'
    sm_Langenau_Mast95@meta.data$treatment = 'NA'
    sm_Langenau_Mast111@meta.data$treatment = 'NA'
    sm_Langenau_Mast139@meta.data$treatment = 'NA'
    sm_Langenau_Mast139_SC@meta.data$treatment = 'NA'
    sm_Langenau_Mast118@meta.data$treatment = 'No'
    SJRHB000026_R2_snRNAseq_small@meta.data$treatment = 'No'
    SJRHB000026_X1_scRNAseq_small@meta.data$treatment = 'No'
    SJRHB000026_R3_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB000026_X2_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB010468_D1_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB010468_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB010927_D1_snRNAseq_small@meta.data$treatment = 'No'
    SJRHB010927_X1_scRNAseq_small@meta.data$treatment = 'No'
    SJRHB010928_R1_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB010928_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB011_D_snRNAseq_small@meta.data$treatment = 'No'
    SJRHB011_X_scRNAseq_small@meta.data$treatment = 'No'
    SJRHB012_R_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB012_Y_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB012_S_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB012_Z_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB012405_D1_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB012405_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013757_D2_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013757_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013758_D1_snRNAseq_small@meta.data$treatment = 'No'
    SJRHB013758_X1_scRNAseq_small@meta.data$treatment = 'No'
    SJRHB013758_D2_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013758_X2_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013759_A1_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013759_X14_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013759_A2_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB013759_X15_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB030680_R1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB030680_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB031320_D1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB031320_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB046156_A1_snRNAseq_small@meta.data$treatment = 'No'
    SJRHB046156_X1_scRNAseq_small@meta.data$treatment = 'No'
    SJRHB049189_D1_snRNAseq_small@meta.data$treatment = 'Yes'
    SJRHB049189_X1_scRNAseq_small@meta.data$treatment = 'Yes'
    Rh41_small@meta.data$treatment = 'NA'
    
    
    
# HOSPITAL  
    PDX104_TS_small@meta.data$hospital = 'Paris'
    PDX29_TS_small@meta.data$hospital = 'Paris'
    PDX35_TS_small@meta.data$hospital = 'Paris'
    Berlin_13304_small@meta.data$hospital = 'Berlin'
    PDX82_small@meta.data$hospital = 'Paris'
    Berlin_13454_small@meta.data$hospital = 'Berlin'
    Berlin_13870_small@meta.data$hospital = 'Berlin'
    Berlin_13933_small@meta.data$hospital = 'Berlin'
    KFR_2_small@meta.data$hospital = 'NA'
    Rh4_2_small@meta.data$hospital = 'NA'
    RMS_2_small@meta.data$hospital = 'NA'
    Rh70_small@meta.data$hospital = 'St. Jude'
    Rh73_small@meta.data$hospital = 'St. Jude'
    Rh71_small@meta.data$hospital = 'St. Jude'
    Rh74_small@meta.data$hospital = 'St. Jude'
    Mast139_small@meta.data$hospital = 'St. Jude'
    Mast118_small@meta.data$hospital = 'St. Jude'
    sm_Langenau_20082@meta.data$hospital = 'MGH'
    sm_Langenau_20696@meta.data$hospital = 'MGH'
    sm_Langenau_21202@meta.data$hospital = 'MGH'
    sm_Langenau_29806@meta.data$hospital = 'MGH'
    sm_Langenau_Mast39@meta.data$hospital = 'St. Jude'
    sm_Langenau_MSK72117@meta.data$hospital = 'MGH'
    sm_Langenau_MSK72117_SC@meta.data$hospital = 'MGH'
    sm_Langenau_MSK74711@meta.data$hospital = 'MGH'
    sm_Langenau_MSK82489@meta.data$hospital = 'MGH'
    sm_Langenau_RD@meta.data$hospital = 'NA'
    sm_Langenau_Mast85_r1@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast85_r2@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast85_r2_SC@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast95@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast111@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast139@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast139_SC@meta.data$hospital = 'St. Jude'
    sm_Langenau_Mast118@meta.data$hospital = 'St. Jude'
    SJRHB000026_R2_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB000026_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB000026_R3_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB000026_X2_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB010468_D1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB010468_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB010927_D1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB010927_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB010928_R1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB010928_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB011_D_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB011_X_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB012_R_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB012_Y_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB012_S_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB012_Z_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB012405_D1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB012405_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013757_D2_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013757_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013758_D1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013758_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013758_D2_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013758_X2_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013759_A1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013759_X14_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013759_A2_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB013759_X15_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB030680_R1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB030680_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB031320_D1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB031320_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB046156_A1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB046156_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB049189_D1_snRNAseq_small@meta.data$hospital = 'St. Jude'
    SJRHB049189_X1_scRNAseq_small@meta.data$hospital = 'St. Jude'
    Rh41_small@meta.data$hospital = 'NA'
     
    
      
##################################################################
############ (5) Combine all datasets ##############
##################################################################
    
    # Saving al objects in RData format
    save(PDX104_TS_small, PDX29_TS_small,  PDX35_TS_small, Berlin_13304_small, 
          PDX82_small, Berlin_13454_small, Berlin_13870_small, Berlin_13933_small, 
          KFR_2_small, Rh4_2_small, RMS_2_small, Rh70_small, 
          Rh73_small,Rh71_small, Rh74_small, 
          Mast139_small,  Mast118_small, sm_Langenau_20082, 
          sm_Langenau_20696, sm_Langenau_21202, sm_Langenau_29806, 
          sm_Langenau_Mast39, sm_Langenau_MSK72117,  sm_Langenau_MSK72117_SC, sm_Langenau_MSK74711, 
          sm_Langenau_MSK82489, sm_Langenau_RD,  sm_Langenau_Mast85_r1,  sm_Langenau_Mast85_r2, sm_Langenau_Mast85_r2_SC, 
          sm_Langenau_Mast95, sm_Langenau_Mast111, sm_Langenau_Mast139,  sm_Langenau_Mast139_SC, 
          sm_Langenau_Mast118,  SJRHB000026_R2_snRNAseq_small,  SJRHB000026_X1_scRNAseq_small, 
          SJRHB000026_R3_snRNAseq_small, SJRHB000026_X2_scRNAseq_small, SJRHB010468_D1_snRNAseq_small, 
          SJRHB010468_X1_scRNAseq_small,  SJRHB010927_D1_snRNAseq_small,  SJRHB010927_X1_scRNAseq_small, 
          SJRHB010928_R1_snRNAseq_small, SJRHB010928_X1_scRNAseq_small, SJRHB011_D_snRNAseq_small, 
          SJRHB011_X_scRNAseq_small,  SJRHB012_R_snRNAseq_small, SJRHB012_Y_scRNAseq_small, 
          SJRHB012_S_snRNAseq_small,  SJRHB012_Z_scRNAseq_small,  SJRHB012405_D1_snRNAseq_small, 
          SJRHB012405_X1_scRNAseq_small,SJRHB013757_D2_snRNAseq_small, SJRHB013757_X1_scRNAseq_small, 
          SJRHB013758_D1_snRNAseq_small,SJRHB013758_X1_scRNAseq_small, SJRHB013758_D2_snRNAseq_small, 
          SJRHB013758_X2_scRNAseq_small,  SJRHB013759_A1_snRNAseq_small, SJRHB013759_X14_scRNAseq_small, 
          SJRHB013759_A2_snRNAseq_small,  SJRHB013759_X15_scRNAseq_small,SJRHB030680_R1_scRNAseq_small,  SJRHB030680_X1_scRNAseq_small, 
          SJRHB031320_D1_scRNAseq_small, SJRHB031320_X1_scRNAseq_small,  SJRHB046156_A1_snRNAseq_small, 
          SJRHB046156_X1_scRNAseq_small,  SJRHB049189_D1_snRNAseq_small, SJRHB049189_X1_scRNAseq_small, Rh41_small,
          file = "/mnt/Sara/write/Danielli_Patel_Langenau_objects.RData")

    ## Combine all samples
    PDX.combined <- merge(PDX104_TS_small, y = list(PDX29_TS_small,  PDX35_TS_small, Berlin_13304_small, 
                                                PDX82_small, Berlin_13454_small, Berlin_13870_small, Berlin_13933_small, 
                                                KFR_2_small, Rh4_2_small, RMS_2_small, Rh70_small, 
                                                Rh73_small,Rh71_small, Rh74_small, 
                                                Mast139_small,  Mast118_small, sm_Langenau_20082, 
                                                sm_Langenau_20696, sm_Langenau_21202, sm_Langenau_29806, 
                                                sm_Langenau_Mast39, sm_Langenau_MSK72117,  sm_Langenau_MSK72117_SC, sm_Langenau_MSK74711, 
                                                sm_Langenau_MSK82489, sm_Langenau_RD,  sm_Langenau_Mast85_r1,  sm_Langenau_Mast85_r2, sm_Langenau_Mast85_r2_SC, 
                                                sm_Langenau_Mast95, sm_Langenau_Mast111, sm_Langenau_Mast139,  sm_Langenau_Mast139_SC, 
                                                sm_Langenau_Mast118,  SJRHB000026_R2_snRNAseq_small,  SJRHB000026_X1_scRNAseq_small, 
                                                SJRHB000026_R3_snRNAseq_small, SJRHB000026_X2_scRNAseq_small, SJRHB010468_D1_snRNAseq_small, 
                                                SJRHB010468_X1_scRNAseq_small,  SJRHB010927_D1_snRNAseq_small,  SJRHB010927_X1_scRNAseq_small, 
                                                SJRHB010928_R1_snRNAseq_small, SJRHB010928_X1_scRNAseq_small, SJRHB011_D_snRNAseq_small, 
                                                SJRHB011_X_scRNAseq_small,  SJRHB012_R_snRNAseq_small, SJRHB012_Y_scRNAseq_small, 
                                                SJRHB012_S_snRNAseq_small,  SJRHB012_Z_scRNAseq_small,  SJRHB012405_D1_snRNAseq_small, 
                                                SJRHB012405_X1_scRNAseq_small,SJRHB013757_D2_snRNAseq_small, SJRHB013757_X1_scRNAseq_small, 
                                                SJRHB013758_D1_snRNAseq_small,SJRHB013758_X1_scRNAseq_small, SJRHB013758_D2_snRNAseq_small, 
                                                SJRHB013758_X2_scRNAseq_small,  SJRHB013759_A1_snRNAseq_small, SJRHB013759_X14_scRNAseq_small, 
                                                SJRHB013759_A2_snRNAseq_small,  SJRHB013759_X15_scRNAseq_small,SJRHB030680_R1_scRNAseq_small,  SJRHB030680_X1_scRNAseq_small, 
                                                SJRHB031320_D1_scRNAseq_small, SJRHB031320_X1_scRNAseq_small,  SJRHB046156_A1_snRNAseq_small, 
                                                SJRHB046156_X1_scRNAseq_small,  SJRHB049189_D1_snRNAseq_small, SJRHB049189_X1_scRNAseq_small, Rh41_small),
                          add.cell.ids = c('PDX104_TS_small',  
                                           'PDX29_TS_small',  
                                           'PDX35_TS_small',  
                                           'Berlin_13304_small',  
                                           'PDX82_small',  
                                           'Berlin_13454_small',  
                                           'Berlin_13870_small',  
                                           'Berlin_13933_small',  
                                           'KFR_2_small',  
                                           'Rh4_2_small',  
                                           'RMS_2_small',  
                                           'Rh70_small',  
                                           'Rh73_small',  
                                           'Rh71_small',  
                                           'Rh74_small',  
                                           'Mast139_small',  
                                           'Mast118_small',  
                                           'sm_Langenau_20082',  
                                           'sm_Langenau_20696',  
                                           'sm_Langenau_21202',  
                                           'sm_Langenau_29806',  
                                           'sm_Langenau_Mast39',  
                                           'sm_Langenau_MSK72117',  
                                           'sm_Langenau_MSK72117_SC',  
                                           'sm_Langenau_MSK74711',  
                                           'sm_Langenau_MSK82489',  
                                           'sm_Langenau_RD',  
                                           'sm_Langenau_Mast85_r1',  
                                           'sm_Langenau_Mast85_r2',  
                                           'sm_Langenau_Mast85_r2_SC',  
                                           'sm_Langenau_Mast95',  
                                           'sm_Langenau_Mast111',  
                                           'sm_Langenau_Mast139',  
                                           'sm_Langenau_Mast139_SC',  
                                           'sm_Langenau_Mast118',  
                                           'SJRHB000026_R2_snRNAseq_small',  
                                           'SJRHB000026_X1_scRNAseq_small',  
                                           'SJRHB000026_R3_snRNAseq_small',  
                                           'SJRHB000026_X2_scRNAseq_small',  
                                           'SJRHB010468_D1_snRNAseq_small',  
                                           'SJRHB010468_X1_scRNAseq_small',  
                                           'SJRHB010927_D1_snRNAseq_small',  
                                           'SJRHB010927_X1_scRNAseq_small',  
                                           'SJRHB010928_R1_snRNAseq_small',  
                                           'SJRHB010928_X1_scRNAseq_small',  
                                           'SJRHB011_D_snRNAseq_small',  
                                           'SJRHB011_X_scRNAseq_small',  
                                           'SJRHB012_R_snRNAseq_small',  
                                           'SJRHB012_Y_scRNAseq_small',  
                                           'SJRHB012_S_snRNAseq_small',  
                                           'SJRHB012_Z_scRNAseq_small',  
                                           'SJRHB012405_D1_snRNAseq_small',  
                                           'SJRHB012405_X1_scRNAseq_small',  
                                           'SJRHB013757_D2_snRNAseq_small',  
                                           'SJRHB013757_X1_scRNAseq_small',  
                                           'SJRHB013758_D1_snRNAseq_small',  
                                           'SJRHB013758_X1_scRNAseq_small',  
                                           'SJRHB013758_D2_snRNAseq_small',  
                                           'SJRHB013758_X2_scRNAseq_small',  
                                           'SJRHB013759_A1_snRNAseq_small',  
                                           'SJRHB013759_X14_scRNAseq_small',  
                                           'SJRHB013759_A2_snRNAseq_small',  
                                           'SJRHB013759_X15_scRNAseq_small',  
                                           'SJRHB030680_R1_scRNAseq_small',  
                                           'SJRHB030680_X1_scRNAseq_small',  
                                           'SJRHB031320_D1_scRNAseq_small',  
                                           'SJRHB031320_X1_scRNAseq_small',  
                                           'SJRHB046156_A1_snRNAseq_small',  
                                           'SJRHB046156_X1_scRNAseq_small',  
                                           'SJRHB049189_D1_snRNAseq_small',  
                                           'SJRHB049189_X1_scRNAseq_small',  
                                           'Rh41_small'),
                          project = "PDX.combo")
    
    ## Reduce dimensions of objects
    DietSeurat(PDX.combined, counts = TRUE, data = TRUE, scale.data = FALSE, features = NULL, assays = 'RNA', dimreducs = FALSE, graphs = FALSE)
    
    
    ### Remove unneeded metadata ###
    PDX.combined[['RNA_snn_res.0.05']] <- NULL
    PDX.combined[['RNA_snn_res.0.3']] <- NULL
    PDX.combined[['RNA_snn_res.0.2']] <- NULL
    PDX.combined[['RNA_snn_res.0.35']] <- NULL
    PDX.combined[['RNA_snn_res.0.4']] <- NULL
    PDX.combined[['RNA_snn_res.0.25']] <- NULL
    PDX.combined[['hash.ID']] <- NULL
    PDX.combined[['ADT_classification']] <- NULL
    PDX.combined[['ADT_classification.global']] <- NULL
    PDX.combined[['ADT_margin']] <- NULL
    PDX.combined[['ADT_secondID']] <- NULL
    PDX.combined[['ADT_maxID']] <- NULL
    PDX.combined[['nFeature_ADT']] <- NULL
    PDX.combined[['nCount_ADT']] <- NULL
    PDX.combined[['lib.size.10k']] <- NULL
    PDX.combined[['RNA_snn_res.0.6']] <- NULL
    PDX.combined[['RNA_snn_res.0.8']] <- NULL
    PDX.combined[['RNA_snn_res.1']] <- NULL
    PDX.combined[['RNA_snn_res.1.2']] <- NULL
    PDX.combined[['RNA_snn_res.1.4']] <- NULL
    PDX.combined[['RNA_snn_res.1.6']] <- NULL
    PDX.combined[['RNA_snn_res.1.8']] <- NULL
    PDX.combined[['RNA_snn_res.2']] <- NULL
    PDX.combined[['counts.human']] <- NULL
    PDX.combined[['fraction.human']] <- NULL
    PDX.combined[['counts.mouse']] <- NULL
    PDX.combined[['fraction.mouse']] <- NULL
    PDX.combined[['sample']] <- NULL
    PDX.combined[['percent.mito']] <- NULL
    PDX.combined[['percent.ribo']] <- NULL
    
    
    # Define mitochondrial genes
    PDX.combined[["percent.mt"]] <- PercentageFeatureSet(PDX.combined, pattern = "^MT-")
    
    
    ###### RMS ###########
    saveRDS(PDX.combined, file = "/mnt/Sara/write/Danielli_Patel_Langenau_20221220.rds")
    

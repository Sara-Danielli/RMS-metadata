# RMS-metadata
Steps used for this project

---

**METADATA ANALYSIS OF ALL RMS SAMPLES (Fig. 1)** 

**1. Merging all RMS datasets from different publications**
- use `0A_merging_RMS.R` script, which:
  - imports .rds files
  - subsets cells to random 1500 cells to speed up computational analysis
  - for Patel et al. the files were mapped to both human and mouse genome; therefore human cells were filtered using 'Tumor' annotation, and only human genes ('hg19') were selected. Gene names were renamed to be consistent with the other publications
  - metadata added
  - datasets combined using 'merge' function

**2. Analyze objects without regressing out inter-sample differences**
- use `1_Integration_no_reg.R` script, to create UMAP plots of samples before regression

**3. Integrate objects using RPCA to remove inter-sample differences**
- use `2_Integration_all_samples_StJude.R` script to remove inter-sample differences (stored in 'name' metacolumn). The script uses RPCA correction (Fig. 1)

**4. Score cells using metaprograms identified in the original publications**
- use `3_Scoring.R` script to score RPCA-integrated object for cell cycle and metaprogram scores; I also tried to score cells using 'UCell' in `3b_Scoring_UCell.R`
- USE `5_Analyses_scores_Vln_plots.R` and `5_Analyses_scores.R` to analyze scoring results (Fig. 3)(if UCell was used: use `4_Analysis_Integration_UCell.R`)

**5. Downstream analysis of RPCA-integrated object**
- For Fig. 1: use `4_Analysis_only.R` script, that:
   - defines high- vs low-cycling cells based on S.scores and G2M.scores
   - runs PCA and UMAP
   - identifies Louvain clusters, their markers, and assigns them names

**6. Visualization plots of RPCA-integrated object**
- Use '4_Plots_only.R' to visualize UMAPS, heatmap, etc.

  ---
 
 **METADATA ANALYSIS OF RMS SAMPLES SPLIT BY SUBTYPE (FN-RMS, FP-RMS PAX3::FOXO1, FP-RMS PAX7::FOXO1) (Fig. 4)** 
 
 **1. Merging all RMS datasets from different publications**
- use `0A_merging_RMS.R` script as explained above
 
**2. Integrate objects using RPCA to remove inter-sample differences and score them**
- use `2_Integration_aRMS_P3F1.R`, '2_Integration_aRMS_P7F1.R' or '2_Integration_eRMS.R' scripts to integrate objects based on molecular subtype for Fig. 4 (FP-RMS (PAX3::FOXO1), FP-RMS (PAX7::FOXO1) or FN-RMS). The script uses RPCA correction.
- the scripts automatically score the objects for module scores

**3. Downstream analysis of RPCA-integrated object**
- For Fig. 4: use `4_Analysis_ARMS_P3F1_only.R`, '4_Analysis_ARMS_P7F1_only.R' or '4_Analysis_ERMS_only.R' scripts to analyze RPCA-integrated objects based on molecular subtype (FP-RMS (PAX3::FOXO1), FP-RMS (PAX7::FOXO1) or FN-RMS). The script:
   - defines high- vs low-cycling cells based on S.scores and G2M.scores
   - runs PCA and UMAP
   - identifies Louvain clusters, their markers, and assigns them names
 
**4. Visualization plots of RPCA-integrated object**
- Use '4_Plots_ARMS_P3F1_only.R', '4_Plots_ARMS_P7F1_only.R' or '4_Plots_ERMS_only.R' to visualize UMAPS, heatmap, etc.
  
**5. Visualize overlap between marker genes across different molecular subgroups (Fig. 3B)**
- Use '10_Correlation_signatures.R' script using the signature markers for each subtype
 
 
     ---
 
 **ANALYSIS OF RMS DUPLICATE SAMPLES SEQUENCED ACROSS DIFFERENT LABS** 

**1. Subset RMS integrated object to samples of interest and plot UMAP and barplot distribution**
- use `7_Analysis_pairs.R` script, which:
  - imports integrated RMS object (from Fig. 1) and subsets it to pairs of interest
  - visualizes distribution of samples on UMAP and bar plot of cell proportions
 
   ---
 
 **ANNOTATION OF RMS TUMOR DATASET BASED ON DEVELOPMENTAL SCRNASEQ DATA** 

**1. Cleanup of Xi et al. Dev. Cell 2020 dataset**
- use `0b_cleanup_Xi_dev_muscle.R` script, which:
  - imports dataframes
  - performs Seurat pipeline for downstream analyses 

**2. Transfer annotations from developmental dataset onto tumors**
- use `8_Labeltransfer_dev_muscle.R` script, which uses SingleR to transfer either developmental time point or annotation onto RMS tumors

   ---
 
 **ANALYSIS OF RMS PATIENT SAMPLE PAIRS OBTAINED DURING DIAGNOSTIC BIOPSY AND DELAYED RESECTION (Fig. 5)** 

**1. FOR ANAND: EXPLAIN HOW RAW DATA WERE ALIGNED**

**2.Score dataset with signature modules and plot scores**
- use `8_Scoring_FFPE_Patient.R` script

 ---
 
 **ANALYSIS OF PDX BIOPSY SAMPLES FROM SJRHB013759_X1 OBTAINED DURING THERAPY (Fig. 5)** 

**1. FOR ANAND: EXPLAIN HOW RAW DATA WERE ALIGNED**

**2.Score dataset with signature modules and plot scores**
- use `9_Scoring_FFPE_treated_aRMS_SJRHB013759_X14.R` script


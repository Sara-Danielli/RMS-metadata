# RMS-metadata
Steps used for this project

---

**METADATA ANALYSIS OF ALL RMS SAMPLES** 

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
- use `2_Integration_all_samples_StJude.R` script to remove inter-sample differences (stored in 'name' metacolumn). The script uses RPCA correction

**4. Score cells using metaprograms identified in the original publications**
- use `3_Scoring.R` script to score RPCA-integrated object for cell cycle and metaprogram scores; I also tried to score cells using 'UCell' in `3b_Scoring_UCell.R`
- USE `5_Analyses_scores_Vln_plots.R` and `5_Analyses_scores.R` to analyze scoring results (if UCell was used: use `4_Analysis_Integration_UCell.R`)

**5. Downstream analysis of RPCA-integrated object**
- use `4_Analysis_Integration.R` script, that:
   - defines cycling cells based on S.scores and G2M. scores > 0
   - defines progenitor and differentiated scores based on mean + 0.25 SD
   - performs Seurat pipeline  
 
  ---
 
 **METADATA ANALYSIS OF RMS SAMPLES SPLIT BY SUBTYPE (FN-RMS, FP-RMS PAX3::FOXO1, FP-RMS PAX7::FOXO1)** 
 
 **1. Merging all RMS datasets from different publications**
- use `0A_merging_RMS.R` script as explained above
 
**2. Integrate objects using RPCA to remove inter-sample differences**
- use `2_Integration_aRMS_P3F1.R` or `2_Integration_aRMS_P7F1.R` or `2_Integration_eRMS.R` to integrate subtype-specific objects and remove inter-sample differences (stored in 'name' metacolumn). The script uses RPCA correction

**3. Downstream analysis of RPCA-integrated object**
- use `4_Analysis_Integration_ARMS_P3F1.R` or `4_Analysis_Integration_ARMS_P7F1.R` or `4_Analysis_Integration_ERMS.R` script which performs same steps es explained above for the integrated RMS object
   
   ---
 
 **ANNOTATION OF RMS TUMOR DATASET BASED ON DEVELOPMENTAL SCRNASEQ DATA** 

**1. Cleanup of Xi et al. Dev. Cell 2020 dataset**
- use `0b_cleanup_Xi_dev_muscle.R` script, which:
  - imports dataframes
  - performs Seurat pipeline for downstream analyses 

**2. Transfer annotations from developmental dataset onto tumors**
- use `8_Labeltransfer_dev_muscle.R` script, which uses SingleR to transfer either developmental time point or annotation onto RMS tumors

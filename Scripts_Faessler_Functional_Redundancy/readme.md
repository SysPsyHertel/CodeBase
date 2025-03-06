# Source codes for 

This repository contains all the necessary files for computing microbiome community models, redundancy measures, and performing statistical analysis on both the  colorectal cancer 
(CRC) and inflammatory bowel disease (IBD) study. 

## Folder Structure

### `CRC_study/`
Contains scripts needed for:

- Computing microbiome community models
- Calculating redundancy measures for the CRC study
- Statistical analysis of the CRC study

### `IBD_study/`
Contains scripts needed for:
- Computing redundancy measures
- Statistical analysis of the IBD studyXxx

## Key Files

### `comparisons_redundancy.R`
Compares the derived measures of functional redundancy with traditional approaches that investigate community-level redundancy (based on Rao's quadratic entropy).

### `compute_redundancy_measures_IBD.R`
Computes redundancy measures using the [FunRed package](https://github.com/SysPsyHertel/FunRed) and prepares the data for statistical analysis in the IBD study.

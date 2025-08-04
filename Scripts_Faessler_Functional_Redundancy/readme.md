### Scripts for "Characterising functional redundancy in microbiome communities via relative entropy"
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.csbj.2025.03.012-blue)](https://doi.org/10.1016/j.csbj.2025.03.012)

This repository contains all the necessary files for computing microbiome community models, redundancy measures, and performing statistical analysis on both the  colorectal cancer 
(CRC) and inflammatory bowel disease (IBD) study. The folder "CRC_study" contains all necessary files for computing microbiome community models, redundancy measures and statistical analysis of the CRCstudy
The folder "IBD_study" contains all necessary files for computing redundancy measures and statistical analysis of the IBD study. The file "comparisons_redundancy.R" compares the derived measures of functional redundancy with traditional appraoches that investigate community-level redundancy (based on Rao's quadratic entropy).

## Files

### `comparisons_redundancy.R`
Compares the derived measures of functional redundancy with traditional approaches that investigate community-level redundancy (based on Rao's quadratic entropy).

### `compute_redundancy_measures_IBD.R`
Computes redundancy measures using the [FunRed package](https://github.com/SysPsyHertel/FunRed) and prepares the data for statistical analysis in the IBD study.


## Folders

### `CRC_study/`
Contains scripts needed for:

- Computing microbiome community models
- Calculating redundancy measures for the CRC study
- Statistical analysis of the CRC study

### `IBD_study/`
Contains scripts needed for:
- Computing redundancy measures
- Statistical analysis of the IBD studyXxx

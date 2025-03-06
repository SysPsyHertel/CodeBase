# Readme
## Files for creating microbiome community models for the colorectal cancer (CRC) study
- run_mgPipe.mat (runs the toolbox)
- JapaneseDiet.txt (input file with diet constraints for creating community models)
- normCoverage_CRC (the normalized abundances from the stem publication, also created in run_mgPipe.mat)
## Files for computing individual secretion fluxes
- compute_individual_secretion_fluxes (microbiome community models needed to be run)
## Files for computing redundancy measures and preparing data for statistical analysis
- metadata.csv (metadata from the stem publication)
- metabolite_data.csv (metabolite data from the stem publication)
- AGORA1_metabolites_abbr.csv (AGORA abbreviations and KEGGIds)
- compute_redundancy_measures_CRC.R (computes redundancy measures using the FunRed package (https://github.com/SysPsyHertel/FunRed) and prepares the data for statistical analysis)
## Files for statistical analysis that includes permutation tests and generation of Supplementary Tables
- statistical_analysis_CRC.R (statistical analysis of the CRC study)

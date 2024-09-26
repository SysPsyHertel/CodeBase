# Readme

## Files for  LP-analyses
- The folder "ComputeLPsolutions" contains files for running the (PSCM) toolbox of the COBRA toolbox. The script "launchIEMAnalysis.m" starts the analyses. It loads the file "Recon3D_S1T8.xlsx", the Supplementary Table S1T8 of Brunk et al. (https://www.nature.com/articles/nbt.4072#Sec25), along with the the file specifying which biomarkers and genes to be tested for knockout, "ST11.xlsx" (see Supplementary Table S11 of the paper).

## Files for curations of whole-body models
- The folder "CurationHarvettaHarvey" contains the files for curation steps to add the metabolite 8-methoxykynurenate and corresponding reactions to Harvey  ("Add_8methoxykynurenate_to_Harvey.m")  and Harvetta ("Add_8methoxykynurenate_to_Harvetta.m").
## Files for creating microbiome community models
- The folder "CreateMicrobiomeModels" contains the necessary files to build microbiome models for individualized host-microbiome models using the normalized abundance table of the CRC dataset of Yachida et al. (https://www.nature.com/articles/s41591-019-0458-7). The script "runMgPipe_CRC.m" starts the pipeline.

## Files for QP-analyses
- The folder "ComputeQPsolutions" contains the files for creating individual microbiome-personalized whole-body models ("createHMmodels.m") and calculating QP solutions of the personalized models, ("getPersonalizedQP_solutions.m") as well as QP solution after knock-out ("getPersonalizedQP_solutions_KO.m"). The script "initializeHMQP.m" starts the pipeline.
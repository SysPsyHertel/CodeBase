CreateMicrobiomeModels contains the necessary files to build microbiome models for individualized host-microbiome models using the normalized abundance table of the CRC dataset of Yachida et al. (https://www.nature.com/articles/s41591-019-0458-7).
The file "runMgPipe_CRC.m" starts the pipeline.

CurationHarvettaHarvey contains the files for curation steps to add the metabolite 8-methoxykynurenate and corresponding reactions to Harvey ("Add_8methoxykynurenate_to_Harvey.m")  and Harvetta ("Add_8methoxykynurenate_to_Harvetta.m").

ComputeQPsolutions contains the files for creating individual microbiome-personalized whole-body models  ("createHMmodels.m") and calculating QP solutions of the personalized models
("getPersonalizedQP_solutions.m") as well as QP solution after knock-out ("getPersonalizedQP_solutions_KO.m"). The file initializeHMQP starts the pipeline.


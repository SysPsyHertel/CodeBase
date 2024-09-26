
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

%% prepare input parameters
% path to and name of the file with abundance information
abunFilePath = [pwd filesep 'normCoverage_CRC.csv'];

% download AGORA
agoraPath = 'C:\Users\MSPG\National University of Ireland, Galway\Group_MSP - Documents\AGORA1\AGORA_1_03\AGORA_1_03_mat';

% number of cores for paralellization
numWorkers = 4;

% create pan-models
panSpeciesModels = [pwd filesep 'panSpeciesModels'];
createPanModels(agoraPath,panSpeciesModels,'Species',numWorkers);


% compute the species to metabolite contributions for metabolites of
% interest% path to file with sample statification information
infoFilePath = [pwd filesep 'metadata.csv'];

% to compute the metabolic profiles (net uptake/secretion) for each sample
computeProfiles = true;

% to ensure models with diet constraints are saved 
saveConstrModels = true;

%% run analysis for the Japanese diet

% path to the file with dietary information
dietFilePath = [pwd filesep 'JapaneseDiet'];

% path to where results will be stored
resPath = [pwd filesep 'MicrobiomeModels_JD' filesep];

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary, statistics, modelsOK] = initMgPipe(panSpeciesModels, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);

metList = {'glutar'};

% define path to models with diet constraints
modPath = [resPath filesep 'Diet'];

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(modPath, 'metList', metList, 'numWorkers', numWorkers);

%% run analysis for the Average European diet

% path to the file with dietary information
dietFilePath='AverageEuropeanDiet';

% path to where results will be stored
resPath=[pwd filesep 'MicrobiomeModels_AD' filesep];

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary, statistics, modelsOK] = initMgPipe(panSpeciesModels, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'saveConstrModels', saveConstrModels, 'numWorkers', numWorkers);

% compute the species to metabolite contributions for metabolites of
% interest
metList = {'glutar'};

% define path to models with diet constraints
modPath = [resPath filesep 'Diet'];

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(modPath, 'metList', metList, 'numWorkers', numWorkers);


abunFilePath = 'C:/Users/faesslerd/Documents/Projects/SHIP_MicrobiomeModelling/MODEL/ag2_mapping_mphl_1200.csv'
[normalizedCoverage,normalizedCoveragePath] = normalizeCoverage(abunFilePath)


initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');

%writecell(normalizedCoverage, 'C:\Users\faesslerd\Documents\Projects\SHIP_MicrobiomeModelling\MODEL\normalizedCoverage.csv')
panSpeciesModels = 'C:\AGORA2\PanModelsAG2'
%createPanModels(agPath, panPath, 'Species')

%agPath = 'C:\AGORA2_01\';
abunFilePath = [pwd filesep 'normalizedCoverage.csv'];
% to compute the metabolic profiles (net uptake/secretion) for each sample
computeProfiles = true;

% path to where results will be stored
resPath = [pwd filesep 'MicrobiomeModels_AED' filesep];

% path to the file with dietary information
dietFilePath = [pwd filesep 'AED'];

% to ensure models with diet constraints are saved 
saveConstrModels = true;

% number of cores for paralellization
numWorkers = 4;

% Start cobratoolbox
[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary, statistics, modelsOK] = initMgPipe(panSpeciesModels, abunFilePath, computeProfiles, 'resPath', resPath, 'dietFilePath', dietFilePath, 'numWorkers', numWorkers);











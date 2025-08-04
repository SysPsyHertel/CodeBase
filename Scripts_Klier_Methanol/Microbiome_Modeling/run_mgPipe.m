modPath = 'C:\AGORA2_01';
abunFilePath = 'C:\Users\faesslerd\Documents\Projects\Methanol_HPM\HMP_models_data\normalizedCoverage.csv';

% path to the file with dietary information
dietFilePath = 'C:\Users\faesslerd\Documents\Projects\Methanol_HPM\mgPIP\AverageEuropeanDiet_4ApplesXyl';

% path to where results will be stored
resPath = 'C:\Users\faesslerd\Documents\Projects\Methanol_HPM\HMP_models_data_4ApplesXyl\';


solverOK=changeCobraSolver('ibm_cplex','LP');
initMgPipe(modPath, abunFilePath, true,  'dietFilePath', dietFilePath,  'resPath', resPath)
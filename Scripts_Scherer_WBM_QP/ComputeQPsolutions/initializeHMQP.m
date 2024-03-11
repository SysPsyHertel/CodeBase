%% Initialize the building of integrated whole-body microbiome models and the QP solutions

% Folder in which the microbiome models are stored
MicromodPath = '/home/openshift/brain/WBM/QP/MicrobiomeModels';

% Path where the integrated whole-body models are saved
HMresPath = '/home/openshift/brain/WBM/QP_RE/New/HM_Models';

% Path where the final tables are stored
resPath = '/home/openshift/brain/WBM/QP_RE/New/Results';

% Load (curated) WBM
load Harvey_104b_8methoxykynurenate_PhysC
male = modelC;

% Create microbiome-personalized WBMs
createHMmodels(male, MicromodPath, HMresPath);

% Get wild-type solutions
getPersonalizedQP_solutions(male, HMresPath,resPath);

% Gene-list for which knock-out QP solutions should be computed
% (VMH Identifier and name)
geneMarkerList = {'5053.1', '', 'PAH';
                   '8942.1', '', 'KYNU'};
               
% Function to perform knockout QP-solutions
getPersonalizedQP_solutions_KO(male, HMresPath, resPath, geneMarkerList);


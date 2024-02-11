% Initialize the building of integrated whole-body microbiome models and the QP solutions

MicromodPath = '/home/openshift/brain/WBM/QP/MicrobiomeModels';
HMresPath = '/home/openshift/brain/WBM/QP_RE/HM_Models';
resPath = '/home/openshift/brain/WBM/QP_RE/104b_Harvey';

% Create microbiome-personalized WBMs
%createHMmodels(MicromodPath, HMresPath);

getPersonalizedQP_solutions(HMresPath,resPath);

geneMarkerList = {'8942.1', '', 'KYNU'};
getPersonalizedQP_solutions_KO(HMresPath, resPath, geneMarkerList);
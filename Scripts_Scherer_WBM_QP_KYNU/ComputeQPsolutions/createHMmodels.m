function  modelHM = createHMmodels(modPath, resPath)
% Function to create host-microbiome WBM models.
%
% INPUT
% modelPath     Path to directory with microbiome models.
% resPath       Path to a directory where host-microbiome WBM models. WBMs are saved
%
% OUTPUT
% Microbiome-personalized Whole-body models
%
% Author:  Tim Hensen, Daniel Fässler, 2023


modelName = 'Harvey';
male = loadPSCMfile(modelName);
model = male;

dInfo = dir(modPath);
modelList={dInfo.name};

for i=3:length(modelList)

microm = load([modPath filesep char(modelList(i))]);

modelF=fieldnames(microm);
microbiota_model=microm.(modelF{1});

% Combine with microbiome model
modelHM = combineHarveyMicrotiota(model, microbiota_model, 400);
modelHM.sex = "female";

% Set model name
modelHM.name=microbiota_model.name;

% Set direction of optimisation
modelHM.osenseStr='min';

% Set microbiome excretion bounds
modelHM.lb(contains(modelHM.rxns,'Excretion_EX_microbiota_LI_biomass'))=1;
modelHM.ub(contains(modelHM.rxns,'Excretion_EX_microbiota_LI_biomass'))=1;

sv = [resPath filesep 'HM_' modelHM.name  '.mat'];
save(sv,'modelHM');
disp(i);
end
end

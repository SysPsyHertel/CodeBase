function  modelHM = createHMmodels(MicromodPath, HMresPath)
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
model = loadPSCMfile(modelName);

dInfo = dir(MicromodPath);
modelList = {dInfo.name};

for i = 3:length(modelList)

    microModel = load(fullfile(MicromodPath, modelList{i}));
    fieldsMicroModel = fieldnames(microModel);
    microbiota_model = microModel.(fieldsMicroModel{1});

    % Combine with microbiome model
    modelHM = combineHarveyMicrotiota(model, microbiota_model, 400);
    modelHM.sex = "male";

    % Set model name and optimization direction
    modelHM.name = microbiota_model.name;
    modelHM.osenseStr = 'min';

    % Set microbiome excretion bounds
    idx = contains(modelHM.rxns, 'Excretion_EX_microbiota_LI_biomass');
    modelHM.lb(idx) = 1;
    modelHM.ub(idx) = 1;

    % Save the model
    sv = fullfile(HMresPath, ['HM_' modelHM.name '.mat']);
    save(sv, 'modelHM');
    disp(i);
end
end
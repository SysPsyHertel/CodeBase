function [qpResults] = getPersonalizedQP_solutions(modPath,resPath)
% Function to optain QP solutions of a personalized WBM.
%
% INPUT
% modelPath     Path to directory with host-microbiome WBM models.
% resPath       Path to a directory where the solution table is saved.
%
% OUTPUT
% Table containing sample names, a variable (feasible) indicating if the
% solution of the QP model exists, solution vector of the QP solution of
% all reactions of the microbiome-personalized WBM.
%
% Author:  Tim Hensen, Daniel Fässler, 2023

changeCobraSolver('ibm_cplex', 'LP', -1);
changeCobraSolver('ibm_cplex', 'QP', -1);

% Obtain all paths to the HM models in the directory
dInfo = dir(modPath);
modelList = string({dInfo(3:end).name});

% Get sample names
sampNames = extractBetween(modelList, 'HM_', '.');

param.minNorm = 1e-6;

% Get rxns from Harvey
modelName = 'Harvey';
male = loadPSCMfile(modelName);

% Preallocate qpResults matrix
qpResults(1,:) = ["sampName"; "feasible"; male.rxns];

for k = 1:length(sampNames)
    % Prepare the variables temporarily storing the simulation results
    model = load(fullfile(modPath, modelList{k}));
    modelF = fieldnames(model);
    model = model.(modelF{1});

    % Fix objective flux bounds at one
    model = changeRxnBounds(model, 'Whole_body_objective_rxn', 1, 'b');

    % Minimize the Euclidean norm of all reactions for a fixed objective
    model.osenseStr = 'min';

    % Set no linear objective
    model.c = zeros(size(model.c));

    % Optimization step
    QP_calc = optimizeWBModel(model, param);

    % Check if solution exists
    if QP_calc.stat == 0
        qpResults(k + 1, :) = [sampNames(k); QP_calc.stat; repelem(NaN, numel(male.rxns))'];
    else
        qpResults(k + 1, :) = [sampNames(k); QP_calc.stat; QP_calc.v(1:numel(male.rxns))];
    end

    disp(k);
end

savepath = [resPath filesep 'WT_Harvey.csv'];
writematrix(qpResults, savepath);
end
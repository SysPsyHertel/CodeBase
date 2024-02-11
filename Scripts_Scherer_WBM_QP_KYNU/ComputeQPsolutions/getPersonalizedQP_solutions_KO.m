function [qpResults_ko] = getPersonalizedQP_solutions_KO_KYNU(modPath,resPath, geneMarkerList)
% Function to optain QP solutions of a personalized WBM after KYNU
% knockout.
%
% INPUT
% modelPath     Path to directory with host-microbiome WBM models.
% resPath       Path to a directory where the solution table is saved.
%
% OUTPUT
% Table containing sample names, a variable (feasible) indicating if the
% solution of the QP model exists, solution vector of the QP solution of
% all reactions of the microbiome-personalized WBM after KYNU knockout.
%
% Author:  Tim Hensen, Daniel Fässler, 2023

causal = 1;

if ~exist('causal', 'var')
    causal = 0;
end

changeCobraSolver('ibm_cplex', 'LP', -1);
changeCobraSolver('ibm_cplex', 'QP', -1);

% Obtain all paths to the HM models in the directory
dInfo = dir(modPath);
modelList = string({dInfo(3:end).name});

% Get sample names
sampNames = extractBetween(modelList, 'HM_', '.');
environment = getEnvironment();

% Set Regularization
param.minNorm = 1e-6;

% Get rxns from Harvey
modelName = 'Harvey';
male = loadPSCMfile(modelName);

% Preallocate qpResults matrix
qpResults_ko(1,:) = ["sampName"; "feasible"; male.rxns];

for g = 1:size(geneMarkerList, 1)
    [IEMRxns, ~] = getRxnsFromGene(male, geneMarkerList{g,1}, causal);

    if ~isempty(IEMRxns)
        for k = 1:length(sampNames)
            restoreEnvironment(environment);
            changeCobraSolver('ibm_cplex', 'LP', 0, -1);
            changeCobraSolver('ibm_cplex', 'QP', 0, -1);

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

            % KO - setting the lower/upper bound of corresponding reactions to zero
            model = changeRxnBounds(model, IEMRxns, 0, 'b');

            % Optimization step
            QP_calc = optimizeWBModel(model, param);

            % Check if solution exists
        if QP_calc.stat == 0
            qpResults_ko(k + 1, :) = [sampNames(k); QP_calc.stat; repelem(NaN, numel(male.rxns))'];
        else
            qpResults_ko(k + 1, :) = [sampNames(k); QP_calc.stat; QP_calc.v(1:numel(male.rxns))];
        end
        end

        savepath = [resPath filesep geneMarkerList{g,3}, '_Harvey.csv'];
        writematrix(qpResults_ko, savepath);
    end
end
end

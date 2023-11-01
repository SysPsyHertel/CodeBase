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


changeCobraSolver('ibm_cplex','LP',-1);
changeCobraSolver('ibm_cplex','QP',-1);

% Obtain all paths to the HM models in the directory
dInfo = dir(modPath);
modelList={dInfo.name};
modelList = string(modelList(1,3:length(modelList)))';

% Get sample names
newA = cellfun(@(x) strsplit(x, 'HM_'), modelList', 'UniformOutput', false);
newA = vertcat(newA{:});
newA = string(newA(:,2));
newA = cellfun(@(x) strsplit(x, '.'), newA, 'UniformOutput', false);
newA = vertcat(newA{:});
sampNames = string(newA(:,1));

param.minNorm = 1e-6;
resultsqp = {};    

if length(sampNames)-1 > 20
    steps=20;
else
    steps=length(sampNames);
end

% Get rxns from Harvey
% In loadPSCMfile of the cobratoolbox needs to be highlighted the actual
% version of the WBM (here Harvey_1_04c)
modelName = 'Harvey';
male = loadPSCMfile(modelName);

qpResults(1,:) = ["sampName"; "feasible"; male.rxns];
lastHumanRxn = length(male.rxns);

for s=1:steps:length(sampNames)
    if length(sampNames)-s>=steps-1
        endPnt=steps-1;
    else
        endPnt=length(sampNames)-s;
    end

    for k=s:s+endPnt
        
    % Prepare the variables temporarily storing the simulation results
    model = load([modPath filesep char(modelList(k))]);
    modelF=fieldnames(model);
    model=model.(modelF{1});
    
	resultsqp{k}{2} = {};
    % Set objective at Whole_body_objective_rxn
    model = changeObjective(model, 'Whole_body_objective_rxn');

    % Fix objective flux bounds at one
    model = changeRxnBounds(model,'Whole_body_objective_rxn',1,'b');
    % Faster with model.c=0;
    % model.c = zeros(size(model.c));
    
    model.sex='male';
    model.osenseStr = 'min';

    % Optimization step
    FBA = optimizeWBModel(model, param);
    % Check if solution exists
	if(FBA.stat == 0)
        qpResultsTmp{k} = [sampNames(k); FBA.stat; repelem(NaN,lastHumanRxn)'];
    else
        qpResultsTmp{k} = [sampNames(k); FBA.stat; FBA.v(1:lastHumanRxn)];
    end
    disp(k);   
    end
    
     for k=s:s+endPnt
        qpResults(k+1,:) = qpResultsTmp{k};
     end
       savepath = [resPath filesep 'NO_KO_Unconstr_Harvey.csv'];
	   writematrix(qpResults, savepath);
end
end

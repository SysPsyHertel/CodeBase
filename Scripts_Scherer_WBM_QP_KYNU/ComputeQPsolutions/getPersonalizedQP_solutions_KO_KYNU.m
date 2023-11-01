function [qpResults_ko] = getPersonalizedQP_solutions_KO_KYNU(modPath,resPath)
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

causal=1
if  ~exist('causal','var')
    causal = 0;
end

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
environment = getEnvironment();

% Set Regularization
param.minNorm = 1e-6;

% Get rxns from Harvey
% In loadPSCMfile of the cobratoolbox needs to be highlighted the actual
% version of the WBM (here Harvey_1_04c)
modelName = 'Harvey';
male = loadPSCMfile(modelName);

qpResults_ko(1,:) = ["sampName"; "feasible"; male.rxns];
lastHumanRxn = length(male.rxns);

% The VMH Identifier of KYNU
geneMarkerList = {'8942.1'}

for g = 1 : size(geneMarkerList,1)
    
    [IEMRxns, grRules] = getRxnsFromGene(male,geneMarkerList{g},causal);
    if(~isempty(IEMRxns))
        for k=1:length(sampNames)
        restoreEnvironment(environment);
        changeCobraSolver('ibm_cplex','LP',0, -1);
        changeCobraSolver('ibm_cplex','QP',0, -1);

        % Prepare the variables temporarily storing the simulation results
        model = load([modPath filesep char(modelList(k))]);
        modelF=fieldnames(model);
        model=model.(modelF{1});
        % Set objective at Whole_body_objective_rxn
        model = changeObjective(model, 'Whole_body_objective_rxn');
        model.sex="male";
        % Fix objective flux bounds at one
        model = changeRxnBounds(model,'Whole_body_objective_rxn',1,'b');

        % Minimise the Euclidean norm of all reactions for a fixed objective
        model.osenseStr = 'min';
        
        %Faster with model.c=0;
        %model.c = zeros(size(model.c));
        
        %KO - setting the lower/upper bound of corresponding reactions to zero
        model = changeRxnBounds(model,IEMRxns,0,'b');   
        disp(k);
        
        % Optimization step
        FBA = optimizeWBModel(model, param);
        
        % Check if solution exists
        if(FBA.stat == 0)
        qpResultsTmp{k} = [sampNames(k); FBA.stat; repelem(NaN,lastHumanRxn)'];
        else
        qpResultsTmp{k} = [sampNames(k); FBA.stat; FBA.v(1:lastHumanRxn)];
        end

        disp("A");   
        end

         for k=1:length(sampNames)
            qpResults_ko(k+1,:) = qpResults_koTmp{k};
         end
         
    savepath = [resPath filesep 'KYNU' 'Harvey_Unconstr.csv'];
	writematrix(qpResults_ko, savepath);
    end
end
end
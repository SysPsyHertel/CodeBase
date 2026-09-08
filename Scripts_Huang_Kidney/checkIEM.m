function [IEMSol] = ...
    checkIEM( ...
    model, IEMRxns, BiomarkerRxns, ...
    minRxnsFluxHealthy, reverseDirObj, ...
    fractionKO, minBiomarker, fixIEMlb, LPSolver)

% =============================================================
% checkIEM
%
% Calculation aligned with the original checkIEM_WBM logic:
%
%   1. Maximize joint WT IEM flux.
%   2. Stop if the joint WT optimum is zero / close to zero.
%   3. Require WT joint flux >= minRxnsFluxHealthy * WT optimum.
%   4. Complete KO of selected IEM reactions.
%   5. Optimize KO joint objective to verify feasibility.
%   6. Maximize each biomarker in constrained WT and KO.
%
% Differences:
%   - optimizeCbModel is used.
%   - no Whole_body_objective_rxn re-optimization step.
%   - WT / KO labels are used.
%
% =============================================================


%% DEFAULTS

if nargin < 4 || isempty(minRxnsFluxHealthy)
    minRxnsFluxHealthy = 1;
end

if nargin < 5 || isempty(reverseDirObj)
    reverseDirObj = 0;
end

if nargin < 6 || isempty(fractionKO)
    fractionKO = 1;
end

if nargin < 7 || isempty(minBiomarker)
    minBiomarker = 0;
end

if nargin < 8 || isempty(fixIEMlb)
    fixIEMlb = 0;
end

if nargin < 9 || isempty(LPSolver)
    LPSolver = 'ibm_cplex';
end


changeCobraSolver(LPSolver,'LP');


tol = 1e-6;

cnt = 1;

IEMSol = {};


%% ============================================================
% CHECK IEM REACTIONS
% ============================================================

if isempty(IEMRxns)
    error('IEMRxns is empty.');
end


for i = 1:length(IEMRxns)

    if ~any(strcmp(model.rxns,IEMRxns{i}))

        error( ...
            'IEM reaction "%s" does not exist in model.', ...
            IEMRxns{i});
    end
end


%% ============================================================
% ADD JOINT IEM OBJECTIVE
% ============================================================

[r,c] = size(model.S);

dummyRxnIdx = c + 1;

dummyMetIdx = r + 1;


model.S(dummyMetIdx,dummyRxnIdx) = -1;


for i = 1:length(IEMRxns)

    idx = ...
        find(strcmp(model.rxns,IEMRxns{i}),1);

    model.S(dummyMetIdx,idx) = 1;
end


model.rxns{dummyRxnIdx,1} = ...
    'IEM_joint_objective_rxn';


model.mets{dummyMetIdx,1} = ...
    'IEM_joint_constraint';


model.lb(dummyRxnIdx,1) = -100000;

model.ub(dummyRxnIdx,1) = 100000;


model.c = zeros(dummyRxnIdx,1);

model.c(dummyRxnIdx) = 1;


if ~isfield(model,'b') || isempty(model.b)

    model.b = zeros(dummyMetIdx,1);

else

    model.b(dummyMetIdx,1) = 0;
end


if ~isfield(model,'csense') || isempty(model.csense)

    model.csense = repmat('E',dummyMetIdx,1);

else

    model.csense(dummyMetIdx,1) = 'E';
end


%% Pad commonly used metadata

if isfield(model,'rxnNames')
    model.rxnNames{dummyRxnIdx,1} = 'IEM joint objective reaction';
end

if isfield(model,'grRules')
    model.grRules{dummyRxnIdx,1} = '';
end

if isfield(model,'rules')
    model.rules{dummyRxnIdx,1} = '';
end

if isfield(model,'rev')
    model.rev(dummyRxnIdx,1) = 1;
end

if isfield(model,'rxnGeneMat')
    model.rxnGeneMat(dummyRxnIdx,:) = ...
        sparse(1,size(model.rxnGeneMat,2));
end

if isfield(model,'metNames')
    model.metNames{dummyMetIdx,1} = 'IEM joint constraint';
end

if isfield(model,'metFormulas')
    model.metFormulas{dummyMetIdx,1} = '';
end

if isfield(model,'metCharges')
    model.metCharges(dummyMetIdx,1) = NaN;
end

if isfield(model,'metComps')

    if isempty(model.metComps)
        model.metComps(dummyMetIdx,1) = 1;
    else
        model.metComps(dummyMetIdx,1) = model.metComps(1);
    end
end


if isfield(model,'C') && ~isempty(model.C)

    model.C = ...
        [model.C, sparse(size(model.C,1),1)];
end


%% ============================================================
% JOINT OBJECTIVE DIRECTION
% ============================================================

if reverseDirObj == 1
    jointDirection = 'min';
else
    jointDirection = 'max';
end


%% ============================================================
% WT JOINT OPTIMUM
% ============================================================

fprintf('\nOptimizing WT joint IEM objective...\n');


solutionWT = ...
    optimizeCbModel(model,jointDirection);


IEMSol{cnt,1} = ...
    'IEM Rxns All obj - WT';


if solutionWT.stat ~= 1

    IEMSol{cnt,2} = 'NA';

    cnt = cnt + 1;


    IEMSol{cnt,1} = ...
        'IEM Rxns All obj - KO';

    IEMSol{cnt,2} = 'NA';

    cnt = cnt + 1;


    for i = 1:size(BiomarkerRxns,1)

        biomarker = BiomarkerRxns{i,1};

        IEMSol{cnt,1} = ['WT:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;

        IEMSol{cnt,1} = ['KO:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;
    end


    return
end


wtIEMFlux = ...
    solutionWT.v(dummyRxnIdx);


if abs(wtIEMFlux) <= tol
    wtIEMFlux = 0;
end


IEMSol{cnt,2} = ...
    num2str(wtIEMFlux);


cnt = cnt + 1;


fprintf( ...
    'WT joint IEM optimum = %.12g\n', ...
    wtIEMFlux);


%% ============================================================
% STOP WHEN JOINT WT FLUX IS ZERO
% ============================================================

if abs(wtIEMFlux) <= tol

    IEMSol{cnt,1} = ...
        'IEM Rxns All obj - KO';

    IEMSol{cnt,2} = 'NA';

    cnt = cnt + 1;


    for i = 1:size(BiomarkerRxns,1)

        biomarker = BiomarkerRxns{i,1};

        IEMSol{cnt,1} = ['WT:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;

        IEMSol{cnt,1} = ['KO:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;
    end


    fprintf('WT joint IEM flux is zero. Biomarker analysis stopped.\n');

    return
end


%% ============================================================
% WT CONSTRAINT
% ============================================================

model.lb(dummyRxnIdx) = ...
    minRxnsFluxHealthy * wtIEMFlux;


model.lb(dummyRxnIdx) = ...
    fix(model.lb(dummyRxnIdx)*1000000)/1000000;


fprintf( ...
    'WT joint lower bound = %.12g (%g x WT optimum)\n', ...
    model.lb(dummyRxnIdx), ...
    minRxnsFluxHealthy);


%% ============================================================
% KO MODEL
% ============================================================

modelKO = model;


if fixIEMlb == 1

    modelKO.lb(dummyRxnIdx) = ...
        (1-fractionKO)*wtIEMFlux;


    modelKO.lb(dummyRxnIdx) = ...
        fix(modelKO.lb(dummyRxnIdx)*1000000)/1000000;

else

    modelKO.lb(dummyRxnIdx) = 0;
end


if fractionKO ~= 1

    error( ...
        ['This implementation follows the complete-KO branch and ' ...
         'requires fractionKO = 1.']);
end


for i = 1:length(IEMRxns)

    idx = ...
        find(strcmp(modelKO.rxns,IEMRxns{i}),1);


    modelKO.lb(idx) = 0;

    modelKO.ub(idx) = 0;
end


%% ============================================================
% KO JOINT OPTIMUM / FEASIBILITY
% ============================================================

fprintf('Optimizing KO joint IEM objective...\n');


solutionKO = ...
    optimizeCbModel(modelKO,jointDirection);


IEMSol{cnt,1} = ...
    'IEM Rxns All obj - KO';


if solutionKO.stat ~= 1

    IEMSol{cnt,2} = 'NA';

    cnt = cnt + 1;


    for i = 1:size(BiomarkerRxns,1)

        biomarker = BiomarkerRxns{i,1};

        IEMSol{cnt,1} = ['WT:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;

        IEMSol{cnt,1} = ['KO:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;
    end


    return
end


koIEMFlux = ...
    solutionKO.v(dummyRxnIdx);


if abs(koIEMFlux) <= tol
    koIEMFlux = 0;
end


IEMSol{cnt,2} = ...
    num2str(koIEMFlux);


cnt = cnt + 1;


fprintf( ...
    'KO joint IEM optimum = %.12g\n', ...
    koIEMFlux);


%% ============================================================
% BIOMARKERS
% ============================================================

for i = 1:size(BiomarkerRxns,1)

    biomarker = ...
        BiomarkerRxns{i,1};


    fprintf('\nTesting biomarker: %s\n',biomarker);


    bioIdx = ...
        find(strcmp(model.rxns,biomarker),1);


    if isempty(bioIdx)

        IEMSol{cnt,1} = ['WT:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;

        IEMSol{cnt,1} = ['KO:' biomarker];
        IEMSol{cnt,2} = 'NA';
        cnt = cnt + 1;

        continue
    end


    %% WT biomarker maximum

    modelWTBio = model;

    modelWTBio.ub(bioIdx) = 100000;

    modelWTBio.c(:) = 0;

    modelWTBio.c(bioIdx) = 1;


    solutionWTBio = ...
        optimizeCbModel(modelWTBio,'max');


    IEMSol{cnt,1} = ...
        ['WT:' biomarker];


    if solutionWTBio.stat == 1

        f = solutionWTBio.v(bioIdx);

        if abs(f) <= tol
            f = 0;
        end

        IEMSol{cnt,2} = num2str(f);

    else

        IEMSol{cnt,2} = 'NA';
    end


    if minBiomarker == 1

        solutionWTMin = ...
            optimizeCbModel(modelWTBio,'min');


        if solutionWTMin.stat == 1

            fmin = solutionWTMin.v(bioIdx);

            if abs(fmin) <= tol
                fmin = 0;
            end

            IEMSol{cnt,4} = num2str(fmin);

        else

            IEMSol{cnt,4} = 'NA';
        end
    end


    cnt = cnt + 1;


    %% KO biomarker maximum

    modelKOBio = modelKO;

    modelKOBio.ub(bioIdx) = 100000;

    modelKOBio.c(:) = 0;

    modelKOBio.c(bioIdx) = 1;


    solutionKOBio = ...
        optimizeCbModel(modelKOBio,'max');


    IEMSol{cnt,1} = ...
        ['KO:' biomarker];


    if solutionKOBio.stat == 1

        f = solutionKOBio.v(bioIdx);

        if abs(f) <= tol
            f = 0;
        end

        IEMSol{cnt,2} = num2str(f);

    else

        IEMSol{cnt,2} = 'NA';
    end


    if minBiomarker == 1

        solutionKOMin = ...
            optimizeCbModel(modelKOBio,'min');


        if solutionKOMin.stat == 1

            fmin = solutionKOMin.v(bioIdx);

            if abs(fmin) <= tol
                fmin = 0;
            end

            IEMSol{cnt,4} = num2str(fmin);

        else

            IEMSol{cnt,4} = 'NA';
        end
    end


    cnt = cnt + 1;
end

end

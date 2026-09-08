function [IEMSolutions, IEMTable, missingMetAll, AccuracySummary] = ...
    performIEMAnalysis_adapted( ...
    model, Recon3D, geneMarkerList, geneMarkerTable, ...
    minRxnsFluxHealthy, causal, reverseDirObj, ...
    fractionKO, minBiomarker, fixIEMlb, LPSolver, ...
    modelLabel)

% =============================================================
% performIEMAnalysis_adapted
%
% IMPORTANT:
%
%   BOTH urine and plasma rows are evaluated using:
%
%       DM_metabolite[c]
%
%   Matrix is NOT used to choose the compartment.
%
%   Matrix is used only for:
%       - keeping urine/plasma observations separate
%       - calculating accuracy separately
%
% Direction:
%
%       KO - WT >  tol -> "+"
%       KO - WT < -tol -> "-"
%       otherwise      -> "="
%
% Existing dnew_complete effect_insilico is preserved as:
%
%       Effect_insilico_Harvey
%
% =============================================================


%% ============================================================
% DEFAULTS
% ============================================================

if nargin < 5 || isempty(minRxnsFluxHealthy)
    minRxnsFluxHealthy = 1;
end

if nargin < 6 || isempty(causal)
    causal = 0;
end

if nargin < 7 || isempty(reverseDirObj)
    reverseDirObj = 0;
end

if nargin < 8 || isempty(fractionKO)
    fractionKO = 1;
end

if nargin < 9 || isempty(minBiomarker)
    minBiomarker = 0;
end

if nargin < 10 || isempty(fixIEMlb)
    fixIEMlb = 0;
end

if nargin < 11 || isempty(LPSolver)
    LPSolver = 'ibm_cplex';
end

if nargin < 12 || isempty(modelLabel)
    modelLabel = 'CurrentModel';
end


directionTol = 1e-6;


%% ============================================================
% CHECK INPUT
% ============================================================

requiredVars = { ...
    'matrix', ...
    'gene', ...
    'biochemical', ...
    'gene_vmhid', ...
    'met_vmhid', ...
    'effect_invivo', ...
    'effect_insilico' ...
    };


if ~istable(geneMarkerTable)
    error('geneMarkerTable must be a MATLAB table.');
end


if ~all(ismember(requiredVars,geneMarkerTable.Properties.VariableNames))
    error('geneMarkerTable is missing required columns.');
end


if size(geneMarkerList,1) ~= height(geneMarkerTable)
    error('geneMarkerList and geneMarkerTable must have the same rows.');
end


%% ============================================================
% SOLVER
% ============================================================

changeCobraSolver(LPSolver,'LP');


%% ============================================================
% FIX BIOMASS MAINTENANCE TO 1 IF PRESENT
% ============================================================

biomassCandidates = { ...
    'biomass_maintenance', ...
    '_biomass_maintenance', ...
    'Whole_body_objective_rxn' ...
    };


biomassIdx = [];


for b = 1:length(biomassCandidates)

    idx = ...
        find(strcmp(model.rxns,biomassCandidates{b}));


    if ~isempty(idx)
        biomassIdx = [biomassIdx; idx(:)];
    end
end


biomassIdx = ...
    unique(biomassIdx);


if ~isempty(biomassIdx)

    model.lb(biomassIdx) = 1;

    model.ub(biomassIdx) = 1;


    fprintf('\nBiomass maintenance fixed to 1:\n');


    for b = 1:length(biomassIdx)

        fprintf( ...
            '  %s\n', ...
            model.rxns{biomassIdx(b)});
    end
end


modelOriginal = ...
    model;


%% ============================================================
% OUTPUT ARRAYS
% ============================================================

nRows = ...
    height(geneMarkerTable);


Model = ...
    repmat(string(modelLabel),nRows,1);


Matrix = ...
    lower(strtrim(string(geneMarkerTable.matrix)));


Gene = ...
    string(geneMarkerTable.gene);


Gene_vmhid = ...
    string(geneMarkerTable.gene_vmhid);


Biochemical = ...
    string(geneMarkerTable.biochemical);


Biomarker = ...
    string(geneMarkerTable.met_vmhid);


Effect_invivo = ...
    strings(nRows,1);


Effect_insilico_Harvey = ...
    strings(nRows,1);


Effect_insilico_CurrentModel = ...
    repmat("NA",nRows,1);


WT = ...
    nan(nRows,1);


KO = ...
    nan(nRows,1);


KO_minus_WT = ...
    nan(nRows,1);


WT_minus_KO = ...
    nan(nRows,1);


Correct_vs_invivo = ...
    nan(nRows,1);


Causal = ...
    repmat(causal,nRows,1);


BiomarkerReaction = ...
    strings(nRows,1);


BiomarkerStatus = ...
    strings(nRows,1);


ReactionStatus = ...
    strings(nRows,1);


AnalysisStatus = ...
    strings(nRows,1);


nRxnsRecon3D = ...
    zeros(nRows,1);


nRxnsModel = ...
    zeros(nRows,1);


nRxnsMissing = ...
    zeros(nRows,1);


ReactionsRecon3D = ...
    strings(nRows,1);


ReactionsModel = ...
    strings(nRows,1);


IEMSolutions = ...
    struct();


missingMetAll = {};


%% ============================================================
% LOOP
% ============================================================

for k = 1:nRows

    model = ...
        modelOriginal;


    Effect_invivo(k) = ...
        normalizeDirection( ...
        geneMarkerTable.effect_invivo(k));


    Effect_insilico_Harvey(k) = ...
        normalizeDirection( ...
        geneMarkerTable.effect_insilico(k));


    geneVMHID = ...
        char(strtrim(Gene_vmhid(k)));


    marker = ...
        char(strtrim(Biomarker(k)));


    matrixType = ...
        Matrix(k);


    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('ROW %d / %d\n',k,nRows);
    fprintf('============================================================\n');

    fprintf('Matrix:     %s\n',matrixType);
    fprintf('Gene:       %s\n',Gene(k));
    fprintf('Gene_vmhid: %s\n',Gene_vmhid(k));
    fprintf('Biomarker:  %s\n',Biomarker(k));
    fprintf('In vivo:    %s\n',Effect_invivo(k));


    %% ========================================================
    % FIND GENE IN RECON3D
    % ========================================================

    reconGene = '';


    idxGene = ...
        find(strcmp(Recon3D.genes,geneVMHID),1);


    if isempty(idxGene)

        geneWithoutDot1 = ...
            regexprep(geneVMHID,'\.1$','');


        idxGene = ...
            find(strcmp( ...
            Recon3D.genes, ...
            geneWithoutDot1), ...
            1);
    end


    if isempty(idxGene) && ...
            ~endsWith(geneVMHID,'.1')

        geneWithDot1 = ...
            [geneVMHID '.1'];


        idxGene = ...
            find(strcmp( ...
            Recon3D.genes, ...
            geneWithDot1), ...
            1);
    end


    if ~isempty(idxGene)

        reconGene = ...
            Recon3D.genes{idxGene};
    end


    %% ========================================================
    % GENE -> REACTION MAPPING
    % ========================================================

    IEMRxnsRecon3D = {};

    IEMRxnsModel = {};

    missingIEMRxns = {};

    reactionMapping = {};

    grRules = {};


    if isempty(reconGene)

        reactionStatus = ...
            'Gene not found in Recon3D';


    else

        [IEMRxnsRecon3D,grRules] = ...
            getRxnsFromGene( ...
            Recon3D, ...
            reconGene, ...
            causal);


        IEMRxnsRecon3D = ...
            unique( ...
            IEMRxnsRecon3D, ...
            'stable');


        for r = 1:length(IEMRxnsRecon3D)

            referenceRxn = ...
                IEMRxnsRecon3D{r};


            [modelRxn,found] = ...
                mapReconReactionToModel( ...
                referenceRxn, ...
                model.rxns);


            if found

                IEMRxnsModel{end+1,1} = ...
                    modelRxn;


                reactionMapping(end+1,:) = ...
                    {referenceRxn,modelRxn};

            else

                missingIEMRxns{end+1,1} = ...
                    referenceRxn;
            end
        end


        IEMRxnsModel = ...
            unique( ...
            IEMRxnsModel, ...
            'stable');


        if isempty(IEMRxnsRecon3D)

            if causal == 1

                reactionStatus = ...
                    'No causal reactions in Recon3D';

            else

                reactionStatus = ...
                    'No associated reactions in Recon3D';
            end


        elseif isempty(IEMRxnsModel)

            if causal == 1

                reactionStatus = ...
                    'No causal reactions in model';

            else

                reactionStatus = ...
                    'No associated reactions in model';
            end


        elseif ~isempty(missingIEMRxns)

            reactionStatus = ...
                'Partial reaction coverage';


        else

            reactionStatus = ...
                'All reactions present';
        end
    end


    ReactionStatus(k) = ...
        string(reactionStatus);


    nRxnsRecon3D(k) = ...
        length(IEMRxnsRecon3D);


    nRxnsModel(k) = ...
        length(IEMRxnsModel);


    nRxnsMissing(k) = ...
        length(missingIEMRxns);


    if isempty(IEMRxnsRecon3D)

        ReactionsRecon3D(k) = "";

    else

        ReactionsRecon3D(k) = ...
            string( ...
            strjoin( ...
            IEMRxnsRecon3D(:)', ...
            ';'));
    end


    if isempty(IEMRxnsModel)

        ReactionsModel(k) = "";

    else

        ReactionsModel(k) = ...
            string( ...
            strjoin( ...
            IEMRxnsModel(:)', ...
            ';'));
    end


    %% ========================================================
    % CYTOSOLIC BIOMARKER FOR BOTH URINE AND PLASMA
    % ========================================================

    metCinput = ...
        [marker '[c]'];


    metMatches = ...
        find(strcmpi( ...
        model.mets, ...
        metCinput));


    BiomarkerRxns = ...
        cell(1,2);


    BiomarkerRxns{1,2} = ...
        'non reported';


    biomarkerAvailable = ...
        false;


    if isempty(metMatches)

        biomarkerStatus = ...
            ['Cytosolic metabolite missing: ' ...
             metCinput];


        missingMetAll{end+1,1} = ...
            metCinput;


        BiomarkerRxns{1,1} = ...
            ['DM_' metCinput];


    elseif length(metMatches) > 1

        error( ...
            'Multiple metabolites matched %s.', ...
            metCinput);


    else

        % Use exact model spelling/capitalization
        metC = ...
            model.mets{metMatches(1)};


        demandRxn = ...
            ['DM_' metC];


        demandIdx = ...
            find(strcmp( ...
            model.rxns, ...
            demandRxn), ...
            1);


        if isempty(demandIdx)

            [model,addedRxns] = ...
                addDemandReaction( ...
                model, ...
                metC, ...
                0);


            demandRxn = ...
                addedRxns{1};


            biomarkerStatus = ...
                'Cytosolic demand reaction added';


        else

            biomarkerStatus = ...
                'Cytosolic demand reaction already present';
        end


        BiomarkerRxns{1,1} = ...
            demandRxn;


        biomarkerAvailable = ...
            true;
    end


    BiomarkerReaction(k) = ...
        string(BiomarkerRxns{1,1});


    BiomarkerStatus(k) = ...
        string(biomarkerStatus);


    fprintf( ...
        'Biomarker reaction: %s\n', ...
        BiomarkerReaction(k));


    fprintf( ...
        'Biomarker status:   %s\n', ...
        BiomarkerStatus(k));


    %% ========================================================
    % RUN IEM
    % ========================================================

    if ~isempty(IEMRxnsModel) && ...
            biomarkerAvailable

        IEMSol = ...
            checkIEM( ...
            model, ...
            IEMRxnsModel, ...
            BiomarkerRxns, ...
            minRxnsFluxHealthy, ...
            reverseDirObj, ...
            fractionKO, ...
            minBiomarker, ...
            fixIEMlb, ...
            LPSolver);


        biomarkerRxn = ...
            BiomarkerRxns{1,1};


        wtRow = ...
            find(strcmp( ...
            IEMSol(:,1), ...
            ['WT:' biomarkerRxn]), ...
            1);


        koRow = ...
            find(strcmp( ...
            IEMSol(:,1), ...
            ['KO:' biomarkerRxn]), ...
            1);


        if ~isempty(wtRow)

            WT(k) = ...
                valueToDouble( ...
                IEMSol{wtRow,2});
        end


        if ~isempty(koRow)

            KO(k) = ...
                valueToDouble( ...
                IEMSol{koRow,2});
        end


        wtJointRow = ...
            find(strcmp( ...
            IEMSol(:,1), ...
            'IEM Rxns All obj - WT'), ...
            1);


        zeroJoint = ...
            false;


        if ~isempty(wtJointRow)

            wtJoint = ...
                valueToDouble( ...
                IEMSol{wtJointRow,2});


            if ~isnan(wtJoint) && ...
                    abs(wtJoint) <= directionTol

                zeroJoint = ...
                    true;
            end
        end


        if zeroJoint

            analysisStatus = ...
                'Not simulated - WT joint IEM flux is zero';


        elseif isnan(WT(k)) || ...
                isnan(KO(k))

            analysisStatus = ...
                'Biomarker optimization unavailable';


        else

            analysisStatus = ...
                'Simulated';
        end


    else

        IEMSol = {};


        if isempty(IEMRxnsModel)

            analysisStatus = ...
                ['Not simulated - ' ...
                 reactionStatus];

        else

            analysisStatus = ...
                ['Not simulated - ' ...
                 biomarkerStatus];
        end
    end


    AnalysisStatus(k) = ...
        string(analysisStatus);


    %% ========================================================
    % DIRECTION
    % ========================================================

    if ~isnan(WT(k)) && ...
            ~isnan(KO(k))

        KO_minus_WT(k) = ...
            KO(k) - WT(k);


        WT_minus_KO(k) = ...
            WT(k) - KO(k);


        Effect_insilico_CurrentModel(k) = ...
            directionFromDelta( ...
            KO_minus_WT(k), ...
            directionTol);


        if isDirection( ...
                Effect_invivo(k))

            Correct_vs_invivo(k) = ...
                double( ...
                Effect_insilico_CurrentModel(k) == ...
                Effect_invivo(k));
        end
    end


    %% ========================================================
    % STORE DETAILS
    % ========================================================

    fieldName = ...
        matlab.lang.makeValidName( ...
        sprintf( ...
        'row_%d_%s_%s_%s', ...
        k, ...
        char(Matrix(k)), ...
        geneVMHID, ...
        marker));


    IEMSolutions.(fieldName).Matrix = ...
        Matrix(k);


    IEMSolutions.(fieldName).Gene = ...
        Gene(k);


    IEMSolutions.(fieldName).Gene_vmhid = ...
        Gene_vmhid(k);


    IEMSolutions.(fieldName).Biochemical = ...
        Biochemical(k);


    IEMSolutions.(fieldName).Biomarker = ...
        Biomarker(k);


    IEMSolutions.(fieldName).BiomarkerReaction = ...
        BiomarkerReaction(k);


    IEMSolutions.(fieldName).Effect_invivo = ...
        Effect_invivo(k);


    IEMSolutions.(fieldName).Effect_insilico_Harvey = ...
        Effect_insilico_Harvey(k);


    IEMSolutions.(fieldName).ReactionsRecon3D = ...
        IEMRxnsRecon3D;


    IEMSolutions.(fieldName).ReactionsModel = ...
        IEMRxnsModel;


    IEMSolutions.(fieldName).ReactionMapping = ...
        reactionMapping;


    IEMSolutions.(fieldName).grRules = ...
        grRules;


    IEMSolutions.(fieldName).solution = ...
        IEMSol;
end


%% ============================================================
% FINAL RESULT TABLE
% ============================================================

IEMTable = ...
    table( ...
    Model, ...
    Matrix, ...
    Gene, ...
    Gene_vmhid, ...
    Biochemical, ...
    Biomarker, ...
    BiomarkerReaction, ...
    Effect_invivo, ...
    Effect_insilico_Harvey, ...
    WT, ...
    KO, ...
    KO_minus_WT, ...
    WT_minus_KO, ...
    Effect_insilico_CurrentModel, ...
    Correct_vs_invivo, ...
    Causal, ...
    ReactionStatus, ...
    BiomarkerStatus, ...
    nRxnsRecon3D, ...
    nRxnsModel, ...
    nRxnsMissing, ...
    ReactionsRecon3D, ...
    ReactionsModel, ...
    AnalysisStatus);


%% ============================================================
% MISSING METABOLITES
% ============================================================

if ~isempty(missingMetAll)

    missingMetAll = ...
        unique( ...
        missingMetAll, ...
        'stable');
end


%% ============================================================
% MATRIX-SPECIFIC ACCURACY
% ============================================================

AccuracySummary = ...
    calculateAccuracy( ...
    IEMTable, ...
    string(modelLabel));


fprintf('\n');
fprintf('============================================================\n');
fprintf('MATRIX-SPECIFIC ACCURACY\n');
fprintf('============================================================\n');


for matrixName = ["urine","plasma"]

    idx = ...
        AccuracySummary.Source == string(modelLabel) & ...
        AccuracySummary.Matrix == matrixName;


    if any(idx)

        x = ...
            AccuracySummary(idx,:);


        fprintf( ...
            '%s - %s: %d / %d correct = %.2f%%\n', ...
            modelLabel, ...
            matrixName, ...
            x.N_correct, ...
            x.N_evaluable, ...
            x.Accuracy_percent);
    end
end

end


% =============================================================
% HELPER: MAP RECON3D REACTION INTO SIMULATION MODEL
% =============================================================
function [modelRxn,found] = ...
    mapReconReactionToModel( ...
    referenceRxn, ...
    modelRxns)

modelRxn = '';

found = ...
    false;


candidates = ...
    {referenceRxn};


if startsWith(referenceRxn,'_')

    if length(referenceRxn) > 1

        candidates{end+1} = ...
            referenceRxn(2:end);
    end


else

    if ~startsWith(referenceRxn,'EX_')

        candidates{end+1} = ...
            ['_' referenceRxn];
    end
end


base = ...
    candidates;


for q = 1:length(base)

    x = ...
        base{q};


    if contains(x,'(e)')

        xE = ...
            strrep( ...
            x, ...
            '(e)', ...
            '[e]');


        candidates{end+1} = ...
            xE;


        candidates{end+1} = ...
            strrep( ...
            xE, ...
            '[e]', ...
            '[d]');
    end


    if contains(x,'[e]')

        candidates{end+1} = ...
            strrep( ...
            x, ...
            '[e]', ...
            '[d]');
    end
end


base = ...
    candidates;


for q = 1:length(base)

    x = ...
        base{q};


    if ~startsWith(x,'_') && ...
            ~startsWith(x,'EX_')

        candidates{end+1} = ...
            ['_' x];
    end
end


candidates = ...
    unique( ...
    candidates, ...
    'stable');


for q = 1:length(candidates)

    idx = ...
        find(strcmp( ...
        modelRxns, ...
        candidates{q}), ...
        1);


    if ~isempty(idx)

        modelRxn = ...
            modelRxns{idx};


        found = ...
            true;


        return
    end
end

end


% =============================================================
% HELPER: CONVERT VALUE TO DOUBLE
% =============================================================
function x = ...
    valueToDouble(value)

if isnumeric(value)

    if isempty(value)

        x = ...
            NaN;

    else

        x = ...
            double(value(1));
    end


    return
end


s = ...
    strtrim(string(value));


if ismissing(s) || ...
        s == "" || ...
        strcmpi(s,"NA") || ...
        strcmpi(s,"NaN") || ...
        strcmpi(s,"ND")

    x = ...
        NaN;


else

    x = ...
        str2double(s);
end

end


% =============================================================
% HELPER: NORMALIZE DIRECTION
% =============================================================
function d = ...
    normalizeDirection(value)

d = ...
    strtrim(string(value));


if ismissing(d) || ...
        d == ""

    d = ...
        "NA";


    return
end


if d == "+" || ...
        d == "-" || ...
        d == "="

    return
end


if strcmpi(d,"NA") || ...
        strcmpi(d,"NaN") || ...
        strcmpi(d,"ND")

    d = ...
        "NA";
end

end


% =============================================================
% HELPER: DELTA -> DIRECTION
% =============================================================
function d = ...
    directionFromDelta( ...
    delta, ...
    tol)

if isnan(delta)

    d = ...
        "NA";


elseif delta > tol

    d = ...
        "+";


elseif delta < -tol

    d = ...
        "-";


else

    d = ...
        "=";
end

end


% =============================================================
% HELPER: VALID DIRECTION?
% =============================================================
function tf = ...
    isDirection(d)

d = ...
    string(d);


tf = ...
    d == "+" || ...
    d == "-" || ...
    d == "=";

end


% =============================================================
% HELPER: ACCURACY
% =============================================================
function Summary = ...
    calculateAccuracy( ...
    IEMTable, ...
    modelLabel)

sources = ...
    [modelLabel; "Harvey"];


matrices = ...
    ["urine";"plasma";"overall"];


rows = {};

cnt = 1;


for s = 1:length(sources)

    source = ...
        sources(s);


    if source == "Harvey"

        prediction = ...
            IEMTable.Effect_insilico_Harvey;


    else

        prediction = ...
            IEMTable.Effect_insilico_CurrentModel;
    end


    for m = 1:length(matrices)

        matrixName = ...
            matrices(m);


        if matrixName == "overall"

            idxMatrix = ...
                true(height(IEMTable),1);


        else

            idxMatrix = ...
                IEMTable.Matrix == matrixName;
        end


        validObserved = ...
            IEMTable.Effect_invivo == "+" | ...
            IEMTable.Effect_invivo == "-" | ...
            IEMTable.Effect_invivo == "=";


        validPredicted = ...
            prediction == "+" | ...
            prediction == "-" | ...
            prediction == "=";


        evaluable = ...
            idxMatrix & ...
            validObserved & ...
            validPredicted;


        N_total = ...
            sum(idxMatrix);


        N_evaluable = ...
            sum(evaluable);


        if N_evaluable > 0

            N_correct = ...
                sum( ...
                prediction(evaluable) == ...
                IEMTable.Effect_invivo(evaluable));


            Accuracy_percent = ...
                100 * ...
                N_correct / ...
                N_evaluable;


        else

            N_correct = ...
                0;


            Accuracy_percent = ...
                NaN;
        end


        rows{cnt,1} = ...
            char(source);


        rows{cnt,2} = ...
            char(matrixName);


        rows{cnt,3} = ...
            N_total;


        rows{cnt,4} = ...
            N_evaluable;


        rows{cnt,5} = ...
            N_correct;


        rows{cnt,6} = ...
            Accuracy_percent;


        cnt = ...
            cnt + 1;
    end
end


Summary = ...
    cell2table( ...
    rows, ...
    'VariableNames',{ ...
    'Source', ...
    'Matrix', ...
    'N_total', ...
    'N_evaluable', ...
    'N_correct', ...
    'Accuracy_percent' ...
    });


Summary.Source = ...
    string(Summary.Source);


Summary.Matrix = ...
    string(Summary.Matrix);

end

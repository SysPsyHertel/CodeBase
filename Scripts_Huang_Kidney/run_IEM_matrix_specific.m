%% ============================================================
% RUN IEM ANALYSIS FOR:
%   1. Bulk PT
%   2. SC PT
%   3. SN PT
%   4. Recon3D
%   5. Male Harvey kidney organ model
%
% Saves ONLY IEMTable for each model
% ============================================================

clc
clear


%% ============================================================
% 1. LOAD dnew_complete.csv
% ============================================================

inputFile = ...
    'C:\Users\faesslerd\Documents\Projects\Chenglong\dnew_complete.csv';

dnew_complete = readtable( ...
    inputFile, ...
    'TextType','string', ...
    'VariableNamingRule','preserve');


%% Keep urine + plasma only

dnew_complete.matrix = ...
    lower(strtrim(string(dnew_complete.matrix)));

dnew_complete = ...
    dnew_complete( ...
        dnew_complete.matrix == "urine" | ...
        dnew_complete.matrix == "plasma", ...
        :);


%% Urine first, plasma second

matrixOrder = ...
    double(dnew_complete.matrix == "plasma");

dnew_complete.MatrixOrder = ...
    matrixOrder;

dnew_complete = ...
    sortrows(dnew_complete,'MatrixOrder');

dnew_complete.MatrixOrder = [];


%% ============================================================
% 2. PREPARE GENE/METABOLITE INPUT
% ============================================================

geneMarkerTable = ...
    dnew_complete;

geneMarkerList = ...
    table2cell( ...
        geneMarkerTable(:,{'gene_vmhid','met_vmhid'}));


%% ============================================================
% 3. LOAD RECON3D REFERENCE
% ============================================================

load Recon3D_Harvey_Used_in_Script_120502

Recon3D_reference = ...
    modelConsistent;


%% ============================================================
% 4. SOLVER / SETTINGS
% ============================================================

changeCobraSolver('ibm_cplex','LP');

minRxnsFluxHealthy = 1;
causal = 1;
reverseDirObj = 0;
fractionKO = 1;
minBiomarker = 0;
fixIEMlb = 0;
LPSolver = 'ibm_cplex';


%% ============================================================
% 5. MODELS
% ============================================================

modelFiles = { ...
    'Bulk_PT_model_constrained_maint_updated_modGPR', ...
    'SC_PT_model_constrained_maint_updated_modGPR', ...
    'SN_PT_model_constrained_modGPR_biomass_maintenance_gapfilled', ...
    'Recon3D', ...
    'OrganAtlas_Harvey' ...
    };


modelLabels = { ...
    'Bulk_PT', ...
    'SC_PT', ...
    'SN_PT', ...
    'Recon3D', ...
    'KidneyOrganMale' ...
    };


%% ============================================================
% 6. OUTPUT DIRECTORY
% ============================================================

outputDir = ...
    'C:\Users\faesslerd\Documents\Projects\Chenglong';


%% ============================================================
% 7. RUN ALL MODELS
% ============================================================

for m = 1:length(modelFiles)

    fprintf('\n\n');
    fprintf('############################################################\n');
    fprintf('MODEL %d / %d: %s\n', ...
        m, ...
        length(modelFiles), ...
        modelLabels{m});
    fprintf('############################################################\n');


    %% ========================================================
    % LOAD SIMULATION MODEL
    % ========================================================

    if strcmp(modelFiles{m},'Recon3D')

        % -----------------------------------------------------
        % Recon3D itself
        % -----------------------------------------------------

        model = ...
            Recon3D_reference;


    elseif strcmp(modelFiles{m},'OrganAtlas_Harvey')

        % -----------------------------------------------------
        % Harvey organ atlas:
        %
        % OrganCompendium_male
        %     .Kidney
        %         .modelAllComp
        % -----------------------------------------------------

        tmp = ...
            load('OrganAtlas_Harvey');


        if ~isfield(tmp,'OrganCompendium_male')

            error( ...
                ['OrganAtlas_Harvey does not contain ' ...
                 'OrganCompendium_male.']);
        end


        if ~isfield(tmp.OrganCompendium_male,'Kidney')

            error( ...
                ['OrganCompendium_male does not contain ' ...
                 'the field Kidney.']);
        end


        if ~isfield( ...
                tmp.OrganCompendium_male.Kidney, ...
                'modelAllComp')

            error( ...
                ['OrganCompendium_male.Kidney does not contain ' ...
                 'modelAllComp.']);
        end


        model = ...
            tmp.OrganCompendium_male.Kidney.modelAllComp;


    else

        % -----------------------------------------------------
        % Kidney models stored in separate MAT files
        % -----------------------------------------------------

        tmp = ...
            load(modelFiles{m});


        % Usually MAT file contains variable "model"
        if isfield(tmp,'model')

            model = ...
                tmp.model;


        else

            % If MAT file contains exactly one variable,
            % use that variable automatically.

            vars = ...
                fieldnames(tmp);


            if length(vars) == 1

                model = ...
                    tmp.(vars{1});


            else

                error( ...
                    ['Could not determine model variable in %s. ' ...
                     'Variables are: %s'], ...
                    modelFiles{m}, ...
                    strjoin(vars,', '));
            end
        end
    end


    %% ========================================================
    % MODEL LABEL
    % ========================================================

    modelLabel = ...
        modelLabels{m};


    %% ========================================================
    % PRINT MODEL INFORMATION
    % ========================================================

    fprintf('\nLoaded model: %s\n',modelLabel);
    fprintf('Reactions:   %d\n',length(model.rxns));
    fprintf('Metabolites: %d\n',length(model.mets));

    if isfield(model,'genes')
        fprintf('Genes:       %d\n',length(model.genes));
    end


    %% ========================================================
    % RUN SAME IEM WORKFLOW
    % ========================================================

    [~, IEMTable, ~, ~] = ...
        performIEMAnalysis_adapted( ...
            model, ...
            Recon3D_reference, ...
            geneMarkerList, ...
            geneMarkerTable, ...
            minRxnsFluxHealthy, ...
            causal, ...
            reverseDirObj, ...
            fractionKO, ...
            minBiomarker, ...
            fixIEMlb, ...
            LPSolver, ...
            modelLabel);


    %% ========================================================
    % SAVE ONLY IEMTable
    % ========================================================

    outputTableFile = ...
        fullfile( ...
            outputDir, ...
            ['IEM_' modelLabel '.xlsx']);


    if exist(outputTableFile,'file')
        delete(outputTableFile);
    end


    writetable( ...
        IEMTable, ...
        outputTableFile);


    fprintf('\nSaved:\n');
    fprintf('%s\n',outputTableFile);

end


%% ============================================================
% DONE
% ============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('ALL MODELS FINISHED\n');
fprintf('============================================================\n');
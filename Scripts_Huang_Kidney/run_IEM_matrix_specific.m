%% ============================================================
% RUN IEM ANALYSIS FOR 3 KIDNEY MODELS + RECON3D
% Saves ONLY IEMTable for each model
% ============================================================

clc
clear


%% ============================================================
% 1. LOAD dnew_complete.csv
% ============================================================

inputFile = ...
    'C:\Users\faesslerd\Documents\Projects\Chenglong\geneMarkerList.csv';

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

Recon3D_reference = modelConsistent;


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
    'Recon3D' ...
    };


modelLabels = { ...
    'Bulk_PT', ...
    'SC_PT', ...
    'SN_PT', ...
    'Recon3D' ...
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


    %% --------------------------------------------------------
    % Load simulation model
    % ---------------------------------------------------------

    if strcmp(modelFiles{m},'Recon3D')

        model = ...
            Recon3D_reference;

    else

        tmp = ...
            load(modelFiles{m});


        % Usually your MAT files contain variable "model"
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


    modelLabel = ...
        modelLabels{m};


    %% --------------------------------------------------------
    % Run exactly the same IEM workflow
    % ---------------------------------------------------------

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


    %% --------------------------------------------------------
    % Save ONLY IEMTable
    % ---------------------------------------------------------

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


    fprintf('\nSaved:\n%s\n',outputTableFile);

end


fprintf('\n');
fprintf('============================================================\n');
fprintf('ALL MODELS FINISHED\n');
fprintf('============================================================\n');
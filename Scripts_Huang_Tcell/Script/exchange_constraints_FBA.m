% ============================================================
% Apply exchange reaction constraints and run FBA on multiple models
% ============================================================


modelFiles = {'A_model.mat', 'C_model.mat', 'B_model.mat'};


modelPath = './models';
outputPath = './output';
constrainFile = './data/constrain.txt';


opts = detectImportOptions(constrainFile, 'Delimiter', '\t');
constrainData = readtable(constrainFile, opts);
rxnIDs = constrainData.Exchange_Reaction_ID;
lbValues = constrainData.Lower_Bound;

for i = 1:length(modelFiles)
    modelFile = modelFiles{i};
    fullModelPath = fullfile(modelPath, modelFile);
    

    loadedModelStruct = load(fullModelPath);
    model = loadedModelStruct.model_out_generic;


    exchangeRxns = findExcRxns(model);
    model.lb(exchangeRxns) = 0;


    for j = 1:length(rxnIDs)
        rxnIndex = find(strcmp(model.rxns, rxnIDs{j}));
        if ~isempty(rxnIndex)
            model.lb(rxnIndex) = lbValues(j);
        end
    end


    model = changeObjective(model, 'biomass_reaction');
    FBAsolution = optimizeCbModel(model, 'max');
    fprintf('Model: %s | Objective Value: %f\n', modelFile, FBAsolution.f);


    newModelName = sprintf('%s_constraint.mat', erase(modelFile, '.mat'));
    save(fullfile(outputPath, newModelName), 'model');
end

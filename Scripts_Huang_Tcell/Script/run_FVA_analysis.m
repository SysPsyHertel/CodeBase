% ============================================================
% Flux Variability Analysis (FVA) for multiple models
% ============================================================


modelFiles = {'A_avg_constrained.mat', 'B_avg_constrained.mat', 'C_avg_constrained.mat'};


modelPath = './models';
outputPath = './output';

for i = 1:length(modelFiles)
    fprintf('Processing model: %s\n', modelFiles{i});
    
  
    loadedData = load(fullfile(modelPath, modelFiles{i}));
    if isfield(loadedData, 'model')
        model = loadedData.model;
    else
        error('Model structure not found in %s', modelFiles{i});
    end
    
 
    model = changeObjective(model, 'biomass_reaction');
    sol = optimizeCbModel(model, 'max');
    if isempty(sol.f) || isnan(sol.f)
        error('No valid biomass_reaction solution for %s', modelFiles{i});
    end
    biomassVal = sol.f;
    fprintf('Optimal biomass for %s: %.6f\n', modelFiles{i}, biomassVal);
    
 
    biomassIdx = find(strcmp(model.rxns, 'biomass_reaction'));
    model.lb(biomassIdx) = 0.95 * biomassVal;
    model.ub(biomassIdx) = biomassVal;
    
  
    solCheck = optimizeCbModel(model, 'max');
    if isempty(solCheck.f) || isnan(solCheck.f)
        error('Model %s is infeasible after biomass constraints. Cannot perform FVA.', modelFiles{i});
    end
    
  
    [minFlux, maxFlux] = fluxVariability(model);
    
  
    totalRxns = length(model.rxns);
    activeRxnIdx = find(abs(minFlux) > 1e-6 | abs(maxFlux) > 1e-6);
    inactiveRxnIdx = find(abs(minFlux) <= 1e-6 & abs(maxFlux) <= 1e-6);
    
    numActive = length(activeRxnIdx);
    numInactive = length(inactiveRxnIdx);
    
    fprintf('Total reactions: %d\n', totalRxns);
    fprintf('Active reactions: %d\n', numActive);
    fprintf('Inactive reactions: %d\n\n', numInactive);
    
  
    T = table(model.rxns, minFlux, maxFlux, ...
        abs(minFlux) > 1e-6 | abs(maxFlux) > 1e-6, ...
        'VariableNames', {'Reaction', 'MinFlux', 'MaxFlux', 'IsActive'});
    
    csvFileName = fullfile(outputPath, [erase(modelFiles{i}, '.mat') '_FVA_results.csv']);
    writetable(T, csvFileName);

end

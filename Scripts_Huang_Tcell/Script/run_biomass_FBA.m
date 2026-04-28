% ============================================================
% Run FBA to compute biomass flux for multiple models
% ============================================================

modelFiles = {'A_avg_constrained.mat', 'B_avg_constrained.mat', 'C_avg_constrained.mat'};

modelPath = '';

for i = 1:length(modelFiles)
  
    modelData = load(fullfile(modelPath, modelFiles{i}));
    model = modelData.model;
    model = changeObjective(model, 'biomass_reaction');
    solution = optimizeCbModel(model, 'max');
    biomassFlux = solution.f;
  
    fprintf('%-40s %-15.6f\n', modelFiles{i}, biomassFlux);
end



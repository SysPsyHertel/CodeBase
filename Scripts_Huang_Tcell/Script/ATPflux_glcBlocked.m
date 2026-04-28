modelFiles = {'A_avg_constrained.mat', 'B_avg_constrained.mat', 'C_avg_constrained.mat'};
modelPath = './Model';

for i = 1:length(modelFiles)
    
    modelData = load([modelPath modelFiles{i}]);
    model = modelData.model;
model = changeRxnBounds(model, 'EX_glc_D[e]', -100, 'l');
        model = changeRxnBounds(model, 'EX_glc_D[e]', 0, 'u'); 

    model = changeObjective(model, 'DM_atp_c_');
    sol1 = optimizeCbModel(model, 'max');
    maxDM_ATP = sol1.f;

    fprintf('Model: %s, max DM_atp_c_ = %.6f\n', modelFiles{i}, maxDM_ATP);

    
    tol = 1e-6;
    rxnID = findRxnIDs(model, 'DM_atp_c_');
    model.lb(rxnID) = maxDM_ATP - tol;  
    model.ub(rxnID) = maxDM_ATP + tol;  


    model = changeObjective(model, 'ATPS4mi');
    sol2 = optimizeCbModel(model, 'max');
    ATPS4mi_val = sol2.x(findRxnIDs(model, 'ATPS4mi'));

    fprintf('Model: %s, ATPS4mi at max DM_atp_c_ = %.6f\n\n', modelFiles{i}, ATPS4mi_val);
end

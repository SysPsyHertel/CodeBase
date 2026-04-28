
modelFiles = {'A_avg_constrained.mat', 'B_avg_constrained.mat', 'C_avg_constrained.mat'};
modelPath = '';


results = {};


results{1,1} = 'Model';
results{1,2} = 'Condition';
results{1,3} = 'Max ATP Flux';
results{1,4} = 'ATPS4mi Flux';

for i = 1:numel(modelFiles)

    data = load(fullfile(modelPath, modelFiles{i}));
    model = data.model;
    

    model = changeObjective(model, 'DM_atp_c_');
    sol = optimizeCbModel(model, 'max');
    results{end+1,1} = modelFiles{i};
    results{end,2}   = 'Baseline (No Constraint)';
    results{end,3}   = sol.f;
    results{end,4}   = sol.x(strcmp(model.rxns,'ATPS4mi'));

    model_O2 = changeRxnBounds(model,'EX_o2[e]',0,'l');
    sol = optimizeCbModel(model_O2, 'max');
    results{end+1,1} = modelFiles{i};
    results{end,2}   = 'O2 uptake blocked';
    results{end,3}   = sol.f;
    results{end,4}   = sol.x(strcmp(model.rxns,'ATPS4mi'));
    

    model = changeRxnBounds(model,'EX_o2[e]',-1000,'l');
    model_lac = changeRxnBounds(model,'EX_lac_L[e]',0,'u');
    sol = optimizeCbModel(model_lac, 'max');
    results{end+1,1} = modelFiles{i};
    results{end,2}   = 'Lactate secretion blocked';
    results{end,3}   = sol.f;
    results{end,4}   = sol.x(strcmp(model.rxns,'ATPS4mi'));
    

    model = changeRxnBounds(model,'EX_lac_L[e]',1000,'u'); 
    model_hex1 = changeRxnBounds(model,'HEX1',0,'b'); 
    sol = optimizeCbModel(model_hex1, 'max');
    results{end+1,1} = modelFiles{i};
    results{end,2}   = 'HEX1 blocked';
    results{end,3}   = sol.f;
    results{end,4}   = sol.x(strcmp(model.rxns,'ATPS4mi'));
end


T = cell2table(results(2:end,:), 'VariableNames', results(1,:));


disp(T)


outFile = fullfile(modelPath, 'Model_ATP_ATPS4mi.xlsx');
writetable(T, outFile);

fprintf('output is saved: %s\n', outFile);

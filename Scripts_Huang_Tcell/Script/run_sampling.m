% ============================================================
% Monte Carlo Sampling for multiple metabolic models
% ============================================================

modelFiles = {'A_avg_constrained.mat', 'B_avg_constrained.mat', 'C_avg_constrained.mat'};


modelPath = './models';
outputPath = './output';
samplingResultPath = './sampling_results';


if ~exist(samplingResultPath, 'dir')
    mkdir(samplingResultPath);
end

for i = 1:length(modelFiles)
  
    modelData = load(fullfile(modelPath, modelFiles{i}));
    model = modelData.model;

    model = changeObjective(model, 'biomass_reaction');

    sol = optimizeCbModel(model, 'max');
    biomassVal = sol.f;
    fprintf('Model: %s | biomass_reaction = %.6f\n', modelFiles{i}, biomassVal);

    biomassIdx = find(strcmp(model.rxns, 'biomass_reaction'));
    model.lb(biomassIdx) = biomassVal;
    model.ub(biomassIdx) = biomassVal;

    warmupn = 10000;               
    fileName = erase(modelFiles{i}, '.mat'); 
    nFiles = 1;                     
    pointsPerFile = 10000;       
    stepsPerPoint = 200;           
    fileBaseNo = 0;                 
    maxTime = 3600000;              
    path = samplingResultPath;      

    fprintf('>>> Start sampling for model: %s\n', modelFiles{i});
    performSampling(model, warmupn, fileName, nFiles, ...
        pointsPerFile, stepsPerPoint, fileBaseNo, maxTime, path);
    fprintf('>>> Sampling finished for model: %s\n\n', modelFiles{i});
end

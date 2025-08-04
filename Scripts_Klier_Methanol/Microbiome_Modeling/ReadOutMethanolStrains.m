modPath = 'C:\AGORA2_01';

% Obtain all paths to the AGORA models in the directory
dInfo = dir(modPath);
modelList={dInfo.name};
modelList = string(modelList(1,3:length(modelList)))';

listspecAG2  = cell(1, 0);
listspecAG2{1,1} = 'MethanolStrainsAG2';

for k=1:length(modelList)
    
    model = load([modPath filesep char(modelList(k))]);
    modelF=fieldnames(model);
    model=model.(modelF{1});

if(~isempty(find(contains(model.rxns, 'PECTIN_DEGe'))))
    listspecAG2{end+1,1} = model.modelName;

end
disp(k);
end

    writecell(listspecAG2, 'Methanol_pect_StrainsAG2.csv')
    
    
    
    
    
    
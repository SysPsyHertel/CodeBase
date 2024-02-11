%% Add 8-Methoxykynurenate to Harvetta_1_04c

%load file corresponding to fileName
modelName = 'Harvetta';

% in loadPSCMfile of the cobratoolbox needs to be highlighted the actual
% version of the WBM (here Harvetta_1_04c)

model = loadPSCMfile(modelName); 

model.rxns(contains(model.rxns, "ahcys"))
model.rxns(contains(model.rxns, "amet"))
model.rxns(contains(model.rxns, "C02470"))

find(contains(model.rxns, "C02470"))
model.mets(find(contains(model.mets, "C02470")))

%%%
% Liver_C02470[c]
% Spleen_C02470[c]
% Adipocytes_C02470[c]
% Retina_C02470[c]
% Heart_C02470[c]
% Ovary_C02470[c]
% Breast_C02470[c]
% Testis_C02470[c]
% Bcells_C02470[c]
% Brain_C02470[c]
% Muscle_C02470[c]
% Agland_C02470[c]
% Uterus_C02470[c]

% Add metabolites in cytosol
modelN = addMetabolite(model, 'Liver_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Spleen_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Adipocytes_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Retina_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Heart_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Ovary_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Breast_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Bcells_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Brain_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Muscle_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Agland_C05830[c]', '8-Methoxykynurenate');
modelN = addMetabolite(modelN, 'Uterus_C05830[c]', '8-Methoxykynurenate');

% Add reaction
model.mets(contains(model.mets, "ahcys"))
model.mets(contains(model.mets, "amet"))

% Adding the reaction containing 8-methoxykynurenate
modelN = addReaction(modelN,'Liver_r03955','reactionFormula','Liver_C02470[c] + Liver_amet[c]  -> Liver_C05830[c] + Liver_ahcys[c]');
modelN = addReaction(modelN,'Spleen_r03955','reactionFormula','Spleen_C02470[c] + Spleen_amet[c]  -> Spleen_C05830[c] + Spleen_ahcys[c]');
modelN = addReaction(modelN,'Adipocytes_r03955','reactionFormula','Adipocytes_C02470[c] + Adipocytes_amet[c]  -> Adipocytes_C05830[c] + Adipocytes_ahcys[c]');
modelN = addReaction(modelN,'Retina_r03955','reactionFormula',' Retina_C02470[c] + Retina_amet[c]  ->  Retina_C05830[c] +  Retina_ahcys[c]');
modelN = addReaction(modelN,'Heart_r03955','reactionFormula','Heart_C02470[c] + Heart_amet[c]  -> Heart_C05830[c] + Heart_ahcys[c]');
modelN = addReaction(modelN,'Ovary_r03955','reactionFormula','Ovary_C02470[c] + Ovary_amet[c]  -> Ovary_C05830[c] + Ovary_ahcys[c]');
modelN = addReaction(modelN,'Breast_r03955','reactionFormula','Breast_C02470[c] + Breast_amet[c]  -> Breast_C05830[c] + Breast_ahcys[c]');
modelN = addReaction(modelN,'Bcells_r03955','reactionFormula','Bcells_C02470[c] + Bcells_amet[c]  -> Bcells_C05830[c] + Bcells_ahcys[c]');
modelN = addReaction(modelN,'Brain_r03955','reactionFormula','Brain_C02470[c] + Brain_amet[c]  -> Brain_C05830[c] + Brain_ahcys[c]');
modelN = addReaction(modelN,'Muscle_r03955','reactionFormula','Muscle_C02470[c] + Muscle_amet[c]  -> Muscle_C05830[c] + Muscle_ahcys[c]');
modelN = addReaction(modelN,'Agland_r03955','reactionFormula','Agland_C02470[c] + Agland_amet[c]  -> Agland_C05830[c] + Agland_ahcys[c]');
modelN = addReaction(modelN,'Uterus_r03955','reactionFormula','Uterus_C02470[c] + Uterus_amet[c]  -> Uterus_C05830[c] + Uterus_ahcys[c]');
modelN = addGenes(modelN,{'8623.1'});
model=modelN;

[n,m] = size(modelN.rxns);
modelN.lb(n-12:n) = -1000000;
modelN.ub(n-12:n) = 1000000;

for i = 1:12
modelN.grRules{n-12+i,1} = '8623.1';
end
model = modelN;

% Add liver exchange
find(string(model.rxns) == 'Liver_r1030')
model.lb(7331)
model.ub(7345)

% Adding metabolites in exchange compartments
model = addMetabolite(model, 'Liver_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Spleen_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Adipocytes_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Retina_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Heart_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Ovary_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Breast_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Bcells_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Brain_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Muscle_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Agland_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Uterus_C05830[e]', '8-Methoxykynurenate');

% Adding exchange reactions
model = addReaction(model,'Liver_tr8mthx','reactionFormula','Liver_C05830[c]  -> Liver_C05830[e]');
model = addReaction(model,'Spleen_tr8mthx','reactionFormula','Spleen_C05830[c]  -> Spleen_C05830[e]');
model = addReaction(model,'Adipocytes_tr8mthx','reactionFormula','Adipocytes_C05830[c]  -> Adipocytes_C05830[e]');
model = addReaction(model,'Retina_tr8mthx','reactionFormula','Retina_C05830[c]  -> Retina_C05830[e]');
model = addReaction(model,'Heart_tr8mthx','reactionFormula','Heart_C05830[c]  -> Heart_C05830[e]');
model = addReaction(model,'Ovary_tr8mthx','reactionFormula','Ovary_C05830[c]  -> Ovary_C05830[e]');
model = addReaction(model,'Breast_tr8mthx','reactionFormula','Breast_C05830[c]  -> Breast_C05830[e]');
model = addReaction(model,'Bcells_tr8mthx','reactionFormula','Bcells_C05830[c]  -> Bcells_C05830[e]');
model = addReaction(model,'Brain_tr8mthx','reactionFormula','Brain_C05830[c]  -> Brain_C05830[e]');
model = addReaction(model,'Muscle_tr8mthx','reactionFormula','Muscle_C05830[c]  -> Muscle_C05830[e]');
model = addReaction(model,'Agland_tr8mthx','reactionFormula','Agland_C05830[c]  -> Agland_C05830[e]');
model = addReaction(model,'Uterus_tr8mthx','reactionFormula','Uterus_C05830[c]  -> Uterus_C05830[e]');

model.lb(n-12:n) = -1000000;
model.ub(n-12:n) = 1000000;
model = addMetabolite(model, 'C05830[bc]', '8-Methoxykynurenate');

% Liver
idx = find(string(model.rxns) == 'Liver_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Liver_EX_C05830(e)_[bc]','reactionFormula','Liver_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Spleen
idx = find(string(model.rxns) == 'Spleen_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Spleen_EX_C05830(e)_[bc]','reactionFormula','Spleen_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Adipocytes
idx = find(string(model.rxns) == 'Adipocytes_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Adipocytes_EX_C05830(e)_[bc]','reactionFormula','Adipocytes_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Retina
idx = find(string(model.rxns) == 'Retina_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Retina_EX_C05830(e)_[bc]','reactionFormula','Retina_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Heart
idx = find(string(model.rxns) == 'Heart_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Heart_EX_C05830(e)_[bc]','reactionFormula','Heart_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Ovary
idx = find(string(model.rxns) == 'Ovary_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Ovary_EX_C05830(e)_[bc]','reactionFormula','Ovary_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Breast
idx = find(string(model.rxns) == 'Breast_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Breast_EX_C05830(e)_[bc]','reactionFormula','Breast_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Bcells
idx = find(string(model.rxns) == 'Bcells_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Bcells_EX_C05830(e)_[bc]','reactionFormula','Bcells_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Muscle
idx = find(string(model.rxns) == 'Muscle_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Muscle_EX_C05830(e)_[bc]','reactionFormula','Muscle_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Agland
idx = find(string(model.rxns) == 'Agland_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Agland_EX_C05830(e)_[bc]','reactionFormula','Agland_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

% Uterus
idx = find(string(model.rxns) == 'Uterus_EX_C02470(e)_[bc]')
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Uterus_EX_C05830(e)_[bc]','reactionFormula','Uterus_C05830[e]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));
% Brain
model = addMetabolite(model, 'C05830[csf]', '8-Methoxykynurenate');

idx = find(string(model.rxns) == 'Brain_EX_C02470(e)_[csf]');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Brain_EX_C05830(e)_[csf]','reactionFormula','Brain_C05830[e]  -> C05830[csf]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

idx = find(string(model.rxns) == 'BBB_C02470[CSF]exp');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'BBB_C05830[CSF]exp','reactionFormula','C05830[csf]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

idx = find(string(model.rxns) == 'BBB_C02470[CSF]upt');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'BBB_C05830[CSF]upt','reactionFormula','C05830[csf]  -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));


%% Kidney
model = addMetabolite(model, 'Kidney_C05830[e]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Kidney_C05830[bcK]', '8-Methoxykynurenate');
model = addMetabolite(model, 'Kidney_C05830[u]', '8-Methoxykynurenate');

idx = find(string(model.rxns) == 'EX_C02470[u]');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'EX_C05830[u]','reactionFormula','C05830[u]  -> ', lowerBound=model.lb(idx), upperBound=model.ub(idx));

idx = find(string(model.rxns) == 'Kidney_EX_C02470(e)_[u]');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Kidney_EX_C05830(e)_[u]','reactionFormula','Kidney_C05830[e] -> C05830[u]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

idx = find(string(model.rxns) == 'Kidney_EX_C02470(e)_[bc]');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Kidney_EX_C05830(e)_[bc]','reactionFormula','Kidney_C05830[e] -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

idx = find(string(model.rxns) == 'Kidney_EX_C02470(e)_[bcK]');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Kidney_EX_C05830(e)_[bcK]','reactionFormula','Kidney_C05830[e] -> Kidney_C05830[bcK]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

idx = find(string(model.rxns) == 'Kidney_EX_C02470[bcK]_[bc]');
model.lb(idx)
model.ub(idx)
model = addReaction(model,'Kidney_EX_C05830[bcK]_[bc]','reactionFormula','Kidney_C05830[bcK] -> C05830[bc]', lowerBound=model.lb(idx), upperBound=model.ub(idx));

save('Harvetta_1_04c_8methoxykynurenate', 'model') 
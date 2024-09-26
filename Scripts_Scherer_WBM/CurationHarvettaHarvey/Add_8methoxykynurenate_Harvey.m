%% Add 8-Methoxykynurenate to Harvey_1_04b

% Set factor for coupling constraints (default 20000)
couplingFactor = 20000;

%load file corresponding to fileName
modelName = 'Harvey';

% in loadPSCMfile of the cobratoolbox needs to be highlighted the actual
% version of the WBM (here Harvey_1_04b)
model = loadPSCMfile(modelName);

% Check which organ contains which reaction
model.rxns(contains(model.rxns, "ahcys"))
model.rxns(contains(model.rxns, "amet"))
model.rxns(contains(model.rxns, "C02470"))

model.mets(find(contains(model.mets, "C02470")))
model.mets(find(contains(model.mets, "amet")))
model.mets(contains(model.mets, "ahcys"))

%%
% Liver_C02470[c]
% Spleen_C02470[c]
% Adipocytes_C02470[c]
% Retina_C02470[c]
% Heart_C02470[c]
% Testis_C02470[c]
% Bcells_C02470[c]
% Brain_C02470[c]
% Muscle_C02470[c]
% Agland_C02470[c]
%%

% Adding the methylation reaction containing 8-methoxykynurenate
model = addReactionsHH(model,{'Liver_r03955'},{'Methylation of 8-methoxykynurenate (Liver)'}, {'Liver_C02470[c] + Liver_amet[c] <=> Liver_C05830[c] + Liver_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Spleen_r03955'},{'Methylation of 8-methoxykynurenate (Spleen)'}, {'Spleen_C02470[c] + Spleen_amet[c] <=> Spleen_C05830[c] + Spleen_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Adipocytes_r03955'},{'Methylation of 8-methoxykynurenate (Adipocytes)'}, {'Adipocytes_C02470[c] + Adipocytes_amet[c] <=> Adipocytes_C05830[c] + Adipocytes_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Retina_r03955'},{'Methylation of 8-methoxykynurenate (Retina)'}, {'Retina_C02470[c] + Retina_amet[c] <=> Retina_C05830[c] +  Retina_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Heart_r03955'},{'Methylation of 8-methoxykynurenate (Heart)'}, {'Heart_C02470[c] + Heart_amet[c] <=> Heart_C05830[c] + Heart_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Testis_r03955'},{'Methylation of 8-methoxykynurenate (Testis)'}, {'Testis_C02470[c] + Testis_amet[c] <=> Testis_C05830[c] + Testis_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Bcells_r03955'},{'Methylation of 8-methoxykynurenate (Bcells)'}, {'Bcells_C02470[c] + Bcells_amet[c] <=> Bcells_C05830[c] + Bcells_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Brain_r03955'},{'Methylation of 8-methoxykynurenate (Brain)'}, {'Brain_C02470[c] + Brain_amet[c] <=> Brain_C05830[c] + Brain_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Muscle_r03955'},{'Methylation of 8-methoxykynurenate (Muscle)'}, {'Muscle_C02470[c] + Muscle_amet[c] <=> Muscle_C05830[c] + Muscle_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);
model = addReactionsHH(model,{'Agland_r03955'},{'Methylation of 8-methoxykynurenate (Agland)'}, {'Agland_C02470[c] + Agland_amet[c] <=> Agland_C05830[c] + Agland_ahcys[c]'},{''},{'Reaction, cytosol'},couplingFactor);

param.SIabsorp = 'n'; % only colon
param.LIabsorp = 'y';% only colon
param.BileDuct = 'n';% only colon
param.Diet = 'n';% only colon
param.PeripheralOrgan = {'Liver', 'Spleen', 'Adipocytes', 'Retina', 'Heart', 'Testis', 'Bcells', 'Brain', 'Muscle', 'Agland'};
Metabolite = 'C05830';
MetaboliteName = '8-methoxykynurenate';
model = addMetFromDiet2Urine2WBM(model,Metabolite,MetaboliteName,param,[],couplingFactor);

sex  = model.sex;
standardPhysiolDefaultParameters;
model = physiologicalConstraintsHMDBbased(model,IndividualParameters);

save('Harvey_104b_8methoxykynurenate_PhysC', 'model')
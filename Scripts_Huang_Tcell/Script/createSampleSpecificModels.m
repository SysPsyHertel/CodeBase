
clear all
setenv('ILOG_CPLEX_PATH','C:\Program Files\IBM\ILOG\CPLEX_Studio129')
initCobraToolbox
solverOK=changeCobraSolver('ibm_cplex','LP');
basePath = '';

normCoverage = readInputTableForPipeline([basePath 'symbol_count_avg.tsv']);


load([basePath 'Mapping_Genes_Recon3D.mat']);



model = readCbModel('Recon3DModel_301.mat');

[C,I] = setdiff(genesMapped(:,3),model.genes,'stable');
genesMapped(I(2:end),:) = [];

rawData = normCoverage(2:end,2:end);
[row, col] = size(rawData);
for i = 1:row
    for j = 1:col
        if ischar(rawData{i,j}) || isstring(rawData{i,j})
            rawData{i,j} = str2double(rawData{i,j});
        end
    end
end
data = cell2mat(rawData);
discretized = discretize_FPKM(data, normCoverage(1,2:end),0);

biomass_rxn = {'biomass_reaction'};


model = changeRxnBounds(model,'biomass_reaction',0.01,'l');

epsilon = 1e-4;
already_mapped_tag = 0;
consensus_proportion = 0.9;


dmem_medium={'aqcobal[e]';'hoxocbl[e]';'ccbl[e]';'C06453[e]';'ala_L[e]';'arg_L[e]';'btn[e]';'ca2[e]';'chol[e]';'cl[e]';'cobalt2[e]';'cu2[e]';'cys_L[e]';'fe2[e]';'fe3[e]';'fol[e]';'glc_D[e]';'gln_L[e]';'gly[e]';'h2o[e]';'hco3[e]';'h2s[e]';'his_L[e]';'ile_L[e]';'inost[e]';'k[e]';'leu_L[e]';'lipoate[e]';'lnlc[e]';'lys_L[e]';'met_L[e]';'mg2[e]';'mn2[e]';'ncam[e]';'o2[e]';'phe_L[e]';'pnto_R[e]';'pro_L[e]';'pydxn[e]';'pyr[e]';'ribflv[e]';'ser_L[e]';'so4[e]';'thf[e]';'thm[e]';'thr_L[e]';'trp_L[e]';'tyr_L[e]';'val_L[e]';'zn2[e]';'co2[e]';'hco3[e]'}; 



unpenalizedSystems = {'Transport, endoplasmic reticular'
    'Transport, extracellular'
    'Transport, golgi apparatus'
    'Transport, mitochondrial'
    'Transport, peroxisomal'
    'Transport, lysosomal'
    'Transport, nuclear'};

subs = {};
for i=1:length(model.subSystems)
    subs{i,1} = model.subSystems{i,1};
end
unpenalized = model.rxns(ismember(subs,unpenalizedSystems));

optional_settings.unpenalized = unpenalized;

optional_settings.func = {'biomass_maintenance','mMMM_cbl','DM_atp_c_','DM_HC00900[m]','MMSAD1m','PPCOACm','r0571','MMALtm','MMMm_cbl','MMEm','MMCDm','METS_FORM','METS_CAT','METS_R_REG','DM_dna5mtc[n]','DM_dna[n]','EX_aqcobal[e]','EX_hoxocbl[e]','EX_ccbl[e]','EX_C06453[e]','EX_ala_L[e]','EX_arg_L[e]','EX_btn[e]','EX_ca2[e]','EX_chol[e]','EX_cl[e]','EX_cobalt2[e]','EX_cu2[e]','EX_cys_L[e]','EX_fe2[e]','EX_fe3[e]','EX_fol[e]','EX_glc_D[e]','EX_gln_L[e]','EX_gly[e]','EX_h2o[e]','EX_hco3[e]','EX_h2s[e]','EX_his_L[e]','EX_ile_L[e]','EX_inost[e]','EX_k[e]','EX_leu_L[e]','EX_lipoate[e]','EX_lnlc[e]','EX_lys_L[e]','EX_met_L[e]','EX_mg2[e]','EX_mn2[e]','EX_ncam[e]','EX_o2[e]','EX_phe_L[e]','EX_pnto_R[e]','EX_pro_L[e]','EX_pydxn[e]','EX_pyr[e]','EX_ribflv[e]','EX_ser_L[e]','EX_so4[e]','EX_thf[e]','EX_thm[e]','EX_thr_L[e]','EX_trp_L[e]','EX_tyr_L[e]','EX_val_L[e]','EX_zn2[e]','EX_co2[e]','EX_hco3[e]'}; % forced additional reactions into the  model
optional_settings.medium = dmem_medium;


mkdir('OutputModels')
Cmodel = model;
samples = normCoverage(1,2:end);


growth = {};
for i = 1:length(samples)
    [model, A_keep] = fastcormics_RNAseq(Cmodel, discretized(:,i), normCoverage(2:end,1), genesMapped, ...
        biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
    model = updateGenes(model);
    

    expressionData = struct;
    expressionData.gene = normCoverage(2:end,1);
    for j=1:length(expressionData.gene)
        findGene = find(strcmp(genesMapped(:,1),expressionData.gene{j}));
        if ~isempty(findGene)
            expressionData.gene{j} = num2str(genesMapped{findGene,3});
        end
    end
    expressionData.value = data(:,i);
    [expressionRxns, parsedGPR, gene_used] = mapExpressionToReactions(model, expressionData);
    expression=struct;
    expression.target=model.rxns;
    expression.value=expressionRxns;
    expression.value(find(isnan(expression.value)))=-1;

    unchRxns = {'biomass_reaction','DM_atp_c_','GTHS','ACONTm', 'AKGDm', 'CSm', 'FUMm', 'ICDHxm', 'ICDHyrm', ...
                          'MDHm', 'SUCD1m', 'SUCOAS1m', 'SUCOASm', 'r0163', 'r0317', ...
                          'r0423', 'r0425', 'r0426', 'r0620', 'ATPS4mi', ...
                          'CYOR_u10mi', 'NADH2_u10mi', 'CYOOm2i', 'FADH2ETC', 'GLYC3PFADm'};
    [C,I] = intersect(expression.target(:,1),unchRxns);
    expression.value(I,1)=-1;
    expression.preprocessed=true;

    model = relaxedApplyEFluxConstraints(model, expression);

    model = changeObjective(model,'biomass_reaction');
    model = changeRxnBounds(model,'DM_dna[n]',0.001,'l');
    model = changeRxnBounds(model,'DM_dna5mtc[n]',0.001,'l');
    
    modelCheck = defineMediumDMEM(model);
    FBA = optimizeCbModel(modelCheck,'max');
    growth{i,1} = samples{i};
    growth{i,2} = FBA.f;

    writeCbModel(model, 'format', 'mat', 'fileName', ['' samples{i} 'TEST.mat']);
end

save growth growth


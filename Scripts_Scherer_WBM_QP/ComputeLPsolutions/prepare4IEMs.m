if ~exist('useSolveCobraLPCPLEX','var')
    global useSolveCobraLPCPLEX
    useSolveCobraLPCPLEX = 0;
end

if ~exist('resultsPath','var')
    global resultsPath
    resultsPath = which('MethodSection3.mlx');
    resultsPath = strrep(resultsPath,'MethodSection3.mlx',['Results' filesep]);
end

changeCobraSolver('ibm_cplex','LP')
changeCobraSolver('ibm_cplex','QP')

if strcmp(modelName,'Harvey')
    %load file corresponding to fileName
    male = loadPSCMfile(modelName);
  
    %standardPhysiolDefaultParameters needs to know what sex it is dealing
    %with
    sex  = male.sex;
    standardPhysiolDefaultParameters;
    
    male = physiologicalConstraintsHMDBbased(male,IndividualParameters);
    EUAverageDietNew;
    %JapaneseDiet;
    male = setDietConstraints(male, Diet);
    model = male;
    modelO = model;
elseif strcmp(modelName,'Harvetta')
    %load file corresponding to fileName
    female = loadPSCMfile(modelName);
    
    %standardPhysiolDefaultParameters needs to know what sex it is dealing
    %with
    sex  = female.sex;
    standardPhysiolDefaultParameters;
    
    female = physiologicalConstraintsHMDBbased(female,IndividualParameters);
    EUAverageDietNew;
    female = setDietConstraints(female, Diet);
    model = female;
    modelO = model;
elseif strcmp(modelName,'Recon3D')
    if useSolveCobraLPCPLEX
        % load Recon3D* and
        load Recon3D_Harvey_Used_in_Script_120502
    else
        modelConsistent = model;
    end
    %makes modifications necessary to adjust Recon3D* for this script
    model = modelConsistent;
    model.rxns = regexprep(model.rxns,'\(e\)','[e]');
    model.rxns = strcat('_',model.rxns);
    model.rxns = regexprep(model.rxns,'_EX_','EX_');
    % add new compartment to Recon
    [model] = createModelNewCompartment(model,'e','u','urine');
    model.rxns = regexprep(model.rxns,'\[e\]_\[u\]','_tr_\[u\]');
    % add exchange reactions for the new [u] metabolites
    U = model.mets(~cellfun(@isempty,strfind(model.mets,'[u]')));
    for i = 1 : length(U)
        model = addExchangeRxn(model,U(i),0,1000);
    end
    
    % create diet reactions
    model.rxns = regexprep(model.rxns,'\[e\]','[d]');
    model.mets = regexprep(model.mets,'\[e\]','[d]');
    EX = model.rxns(~cellfun(@isempty,strfind(model.rxns,'EX_')));
    D = model.rxns(~cellfun(@isempty,strfind(model.rxns,'[d]')));
    EX_D = intersect(EX,D);
    model.rxns(ismember(model.rxns,EX_D)) = strcat('Diet_', model.rxns(ismember(model.rxns,EX_D)));
    % apply diet constraints
    model = setDietConstraints(model);
    % only force lower constraints of diet as no fecal outlet
    EX_D = model.rxns(strmatch('Diet_',model.rxns));
    model.ub(ismember(model.rxns,EX_D)) = 0;
    model.lb(ismember(model.rxns,'Diet_EX_o2[d]')) = -1000;
    
    if useSolveCobraLPCPLEX
        model.A = model.S;
    else
        if isfield(model,'A')
            model = rmfield(model,'A');
        end
    end
    model.rxns(ismember(model.rxns,'_biomass_maintenance')) = {'Whole_body_objective_rxn'};
    model.lb(ismember(model.rxns,'Whole_body_objective_rxn')) = 1;
    model.ub(ismember(model.rxns,'Whole_body_objective_rxn')) = 1;
    modelO = model;
end
cnt = 1;
minRxnsFluxHealthy = 1;%0.9;

%% integrate microbes into the whole-body reconstructions
% load microbe model
microbiome = 0;
set = 0;
if microbiome == 1
    files = {'SRS011239'};
    if set == 1
        S= load(strcat('/microbiota_model_samp_',files{1},'.mat'));
    end
    microbiota_model = S.microbiota_model;
    microbiota_model.rxns = strcat('Micro_',microbiota_model.rxns);
    modelHM = combineHarveyMicrotiota(model,microbiota_model,400);
    
    modelHM = changeRxnBounds(modelHM,'Whole_body_objective_rxn',1,'b');
    modelHMO = modelHM;
    % now set all strains to 0 but 1
    %Bacteroides_thetaiotaomicron_VPI_5482_biomass[c]
    % modelHM.S = modelHM.A;
    % modelHM.S(:,strmatch('communityBiomass',modelHM.rxns))=0;
    % modelHM.S(strmatch('Bacteroides_thetaiotaomicron_VPI_5482_biomass[c]',modelHM.mets),strmatch('communityBiomass',modelHM.rxns))=-1;
    
    % modelHM.S(strmatch('microbiota_LI_biomass[luM]',modelHM.mets),strmatch('communityBiomass',modelHM.rxns))=1;
    % modelHM.A = modelHM.S;
    
    modelHM.lb(ismember(modelHM.rxns,'Excretion_EX_microbiota_LI_biomass[fe]'))=0.1; %
    modelHM.ub(ismember(modelHM.rxns,'Excretion_EX_microbiota_LI_biomass[fe]'))=1; %
    
    model = modelHM;
end

%% set unified reaction constraints -- they are duplicated again in individual scripts

R = {'_ARGSL';'_GACMTRc';'_FUM';'_FUMm';'_HMR_7698';'_UAG4E';'_UDPG4E';'_GALT'; '_G6PDH2c';'_G6PDH2r';'_G6PDH2rer';...
    '_GLUTCOADHm';'_r0541'; '_ACOAD8m';'_RE2410C';'_RE2410N'};
RxnsAll2 = '';
for i = 1: length(R)
    RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R{i})));
    RxnsAll2 =[RxnsAll2;RxnsAll];
end

%excluded reactions
R2 = {'_FUMt';'_FUMAC';'_FUMS';'BBB'};
RxnsAll4 = '';
for i = 1: length(R2)
    RxnsAll3 = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
    RxnsAll4 =[RxnsAll4;RxnsAll3];
end
RxnsAll4 = unique(RxnsAll4);
IEMRxns = setdiff(RxnsAll2,RxnsAll4);
RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
if ~isempty(RxnMic)
    RxnMic
end
IEMRxns = setdiff(IEMRxns,RxnMic);
% set ARGSL to be irreversible
model.lb(ismember(model.rxns,IEMRxns)) = 0;

R2 = {'_r0784';'_r0463'};
RxnsAll2 = '';
for i = 1: length(R2)
    RxnsAll = model.rxns(~cellfun(@isempty,strfind(model.rxns,R2{i})));
    RxnsAll2 =[RxnsAll2;RxnsAll];
end
X = unique(RxnsAll2);
RxnMic = model.rxns(~cellfun(@isempty,strfind(model.rxns,'Micro_')));
if ~isempty(RxnMic)
    RxnMic
end
X = setdiff(X,RxnMic);
model.lb(ismember(model.rxns,X)) = 0;
model.ub(ismember(model.rxns,X)) = 0;

%%%
Rnew = {'BileDuct_EX_12dhchol[bd]_[luSI]';'BileDuct_EX_3dhcdchol[bd]_[luSI]';'BileDuct_EX_3dhchol[bd]_[luSI]';'BileDuct_EX_3dhdchol[bd]_[luSI]';'BileDuct_EX_3dhlchol[bd]_[luSI]';'BileDuct_EX_7dhcdchol[bd]_[luSI]';'BileDuct_EX_7dhchol[bd]_[luSI]';'BileDuct_EX_cdca24g[bd]_[luSI]';'BileDuct_EX_cdca3g[bd]_[luSI]';'BileDuct_EX_cholate[bd]_[luSI]';'BileDuct_EX_dca24g[bd]_[luSI]';'BileDuct_EX_dca3g[bd]_[luSI]';'BileDuct_EX_dchac[bd]_[luSI]';'BileDuct_EX_dgchol[bd]_[luSI]';'BileDuct_EX_gchola[bd]_[luSI]';'BileDuct_EX_hca24g[bd]_[luSI]';'BileDuct_EX_hca6g[bd]_[luSI]';'BileDuct_EX_hdca24g[bd]_[luSI]';'BileDuct_EX_hdca6g[bd]_[luSI]';'BileDuct_EX_hyochol[bd]_[luSI]';'BileDuct_EX_icdchol[bd]_[luSI]';'BileDuct_EX_isochol[bd]_[luSI]';'BileDuct_EX_lca24g[bd]_[luSI]';'BileDuct_EX_tchola[bd]_[luSI]';'BileDuct_EX_tdchola[bd]_[luSI]';'BileDuct_EX_tdechola[bd]_[luSI]';'BileDuct_EX_thyochol[bd]_[luSI]';'BileDuct_EX_uchol[bd]_[luSI]'};
model.ub(ismember(model.rxns,Rnew)) = 100;

modelO = model;
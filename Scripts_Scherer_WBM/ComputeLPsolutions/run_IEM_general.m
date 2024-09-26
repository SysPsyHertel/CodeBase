% This script predicts known biomarker metabolites in
% different biofluid compartments (urine, blood, csf) of the whole-body
% model for 57 inborn-errors of metabolism (IEMs).
% The meaning of the abbreviations for metabolites and IEMs used in this
% script can be found at www.vmh.life.
% The supported model options are 'male', 'female', and 'Recon3D'. Please
% define those using the variable 'sex' (e.g., sex = 'male').
%
% Ines Thiele 2018 - 2019

modelName = 'Harvey';
%modelName = 'Harvetta';
urine = 1; % test for urine metabolites
compartment = {'[bc]'};
%% original gene and biomarker list (Fall 2021)
if 0
    geneMarkerList = {
        %   %  gene  biomarker list
        '249.1' '3pg;cholp;glyc3p;ethamp'
        '1491.1' 'cyst_L' %??
        '56267.1' 'im4ac'
        '1806.1'	'thym;ura'
        '51179.1' 'r34hpl;indlac'
        '8942.1'	'C02470'
        '95.1'	'CE1554;CE1556;acglu;acgly;acile_L;acleu_L;C02712;acser;acthr_L'
        '56954.1' 'HC00591'
        '340024.1' 'asn_L;his_L'
        '29958.1'	'dmgly'
        '6584.1'	'crn'
        '435.1'	'argsuc'
        '6561.1'	'so4'
        '65263.1'	'1pyr5c'
        '26873.1'	'5oxpro'
        '883.1'	'r34hpl;indlac;phpyr' % r34hpl is microbial origin
        '1757.1'	'sarcs'
        %    '1757.1'	   'c4dc;HC00900;hmcarn;HC00342;succ;4abut;urate;ppbng;dchac;dca24g;dca3g;co;pheme;C05767;C05770;ppp9'
        '10249.1'	'tiggly'
        '6539.1'	'glyb'
        '3034.1'	'his_L;urcan'
        '5053.1'	'phe_L;phlac;phpyr'
        '35.1'	'ethmalac'
        '4329.1'	'3hmp'
        '6652.1'	'rbt'
        '123876.1'	'C10164'
        '197322.1'	'ethmalac;HC00900'
        '125061.1'	'nformanth'
        '84735.1'	'carn;hmcarn'
        '11136.1'	'pcollg5hlys;cys_L;lys_L;orn'
        '10841.1'	'forglu'
        '51733.1'	'3uib;cala'
        %
        '51733.1'  '56dura;56dthm;cala;3uib'
        '3155.1' '3hivac;3mglutac'
        
        '2629.1' 'gluside_hs'
        };
end

%%
if 0
    % load WBM models, make adjustment to WBM models, as in publication, and also define the
    % required inputs to performIEMAnalysis
    prepare4IEMs;
    model = generateWBM_104(model);
    modelori = model;
else
    if strcmp(modelName,'Harvey')
        load('Harvey_104b.mat')
    elseif strcmp(modelName,'Harvetta')
        load('Harvetta_1_04c.mat')
    end
    prepare4IEMs;
    modelori = model;
end
%% these changes where not part of the first Johannes Freiburg paper

model = modelori;
model.ub(contains(model.rxns,'HMR_4701'))=0;
model.lb(contains(model.rxns,'HMR_4701'))=0;
model.lb(contains(model.rxns,'GLYSARCNc'))=0;
model.ub(contains(model.rxns,'HMR_4700'))=0;
model.lb(contains(model.rxns,'HMR_4700'))=0;
% reaction DURAD and DURAD2 should go in reverse direction
model.S(:,contains(model.rxns,'DURAD'))=-1*model.S(:,contains(model.rxns,'DURAD'));
model.ub(contains(model.rxns,'DURAD'))=100000;
model.lb(contains(model.rxns,'DURAD'))=0;
model.ub(contains(model.rxns,'r0236'))=0;
model.lb(contains(model.rxns,'r0236'))=0;
model.lb(contains(model.rxns,'CHAT'))=0;% HMDB: choline is essential vit and precursor of ach
model.lb(contains(model.rxns,'PSSA1_hs'))=-1; % HMDB: only small reverse flux

%model.ub(contains(model.rxns,'HMR_9604'))=0;
%model.lb(contains(model.rxns,'HMR_9604'))=0;
% this is not need anymore with the change in ko'ing the objective however
% the issue of inconsistent transporters remains
% model.ub(contains(model.rxns,'GLYBt4_2_r'))=0;
%model.lb(contains(model.rxns,'GLYBt4_2_r'))=0; % this should be the correct transporter however, it cannot transport glyb anymore if I delete it
% model.lb(contains(model.rxns,'BETBGTtc'))=0; % hence I keep this one for the moment but it should be corrected (TODO)
% model.ub(contains(model.rxns,'BETBGTtc'))=0;
% model.lb(contains(model.rxns,'GABABGTtc'))=0;
% model.ub(contains(model.rxns,'GABABGTtc'))=0;

% remove inconsistent reaction from gene 6584.1
% unfortunately, CRNtuNa is again the correct reaction but most organs have
% CRNt - needs to be fixes (TODO)
causal = 1;
[Rxns, grRules] = getRxnsFromGene(model,'6584.1',causal);
Rxns=Rxns(contains(Rxns,'_CRNtuNa'));
model.lb(contains(model.rxns,Rxns))=0;
model.ub(contains(model.rxns,Rxns))=0;

model.ub(contains(model.rxns,'HMR_4700'))=0;

% add indole to diet
model.lb(contains(model.rxns,'Diet_EX_indole[d]'))=-1;

%% This runs the original version that I sent to Johannes and Freiburg
if 1
    clear IEMSolutions IEMTable missingMet
    % [IEMSolutions,IEMTable,missingMet] = performIEMAnalysis(model,geneMarkerList,compartment,urine,minRxnsFluxHealthy, reverseDirObj, fractionKO,minBiomarker,fixIEMlb, LPSolver)
    % anything that is not defined will be used as in defaults
    causal = 1;
    minRxnsFluxHealthy=1
    [IEMSolutions_causal,IEMTable2_causal,missingMet2] = performIEMAnalysis(model,geneMarkerList,compartment,urine,minRxnsFluxHealthy,causal);
    %causal = 0
    %[IEMSolutions_noncausal,IEMTable2_noncausal,missingMet2] = performIEMAnalysis(model,geneMarkerListMarkerList,compartment,urine,minRxnsFluxHealthy,causal);
end
%% here I am investigating the use of Duals for finding known and novel biomarkers
% It appears that the original version (above) can pick up biomarkers that
% are not limited by governing constraints, this method allows for finding
% some of the known and novel ones (as long as they are constraint by governing
% constraints) - IT Jan 2022

clear IEMBiomarkers*
if 0
    for k = 1 : size(geneMarkerList,1)
        [IEMRxns, grRules] = getRxnsFromGene(model,geneMarkerList{k},causal);
        [ PotBioMarkersIEM ,PotBioMarkersIEMR,maxSol] = getPotMarkersIEM(model,IEMRxns,'new');
        G = ['G_' geneMarkerList{k}];
        G = regexprep(G,'\.','_');
        IEMBiomarkers.(G).PotBioMarkersIEM = PotBioMarkersIEM;
        IEMBiomarkers.(G).PotBioMarkersIEMR = PotBioMarkersIEMR;
        IEMBiomarkers.(G).IEMRxns = IEMRxns;
        IEMBiomarkers.(G).maxSol = maxSol;
        IEMBiomarkers.(G).IEMRxnsgrRules = grRules;
        IEMBiomarkers.(G).IEMRxnsgrRulesUnique = unique(grRules);
    end
    
    if 0 %- this method did not work too good
        % I seem to loose hits when I am doing this
        bc_mets = model.mets(contains(model.mets,'[bc]'));
        model_exp = model;
        % add demand reaction to each bc metabolite
        [model_exp] = addDemandReaction(model_exp,bc_mets,0);
        
        for k = 1 : size(geneMarkerList,1)
            [IEMRxns, grRules] = getRxnsFromGene(model_exp,geneMarkerList{k,1},causal);
            [ PotBioMarkersIEM ,PotBioMarkersIEMR]= getPotMarkersIEM(model_exp,IEMRxns);
            G = ['G_' geneMarkerList{k,1}];
            G = regexprep(G,'\.','_');
            IEMBiomarkers_exp.(G).PotBioMarkersIEM = PotBioMarkersIEM;
            IEMBiomarkers_exp.(G).PotBioMarkersIEMR = PotBioMarkersIEMR;
            IEMBiomarkers_exp.(G).IEMRxns = IEMRxns;
            IEMBiomarkers_exp.(G).IEMRxnsgrRulesUnique = unique(grRules);
            IEMBiomarkers_exp.(G).IEMRxnsgrRules = grRules;
        end
    end
end
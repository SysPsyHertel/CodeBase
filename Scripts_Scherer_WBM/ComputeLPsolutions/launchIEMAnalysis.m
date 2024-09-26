%read in Freiburg data
if 1
    filename = 'ST11.xlsx';
else
    filename = 'table';
    filename = 'Supplementary Tables_submitT3.xlsx';
end
[NUM,TXT,RAW]=xlsread(filename);
InputTable = RAW;

% Load Recon 3D gene information
filenameRecon3D = 'Recon3D_S1T8.xlsx';
%[NUM,TXT,RAW]=xlsread(filenameRecon3D,'Supplement Table S9');
[NUM,TXT,RAW]=xlsread(filenameRecon3D);%,'Supplement Table S9');
Recon3DGeneInfo = RAW(1:end,1:19);

% find EntrezGene ID column in Recon3DGeneInfo
colE = find(contains(Recon3DGeneInfo(1,:),'EntrezGene'));
% find HGNC ID column in Recon3DGeneInfo
colHGNC = find(contains(Recon3DGeneInfo(1,:),'HGNC'));

% add missing mappings to Recon3DGeneInfo - due to different gene symbols
% being used
[r,c] = size(Recon3DGeneInfo);
missingGenes = {'KYAT3' '56267'
    'SLC17A1/A3/A4'  '6568' % not clear which gene it is I assume SLC17A1 here
    };
for i = 1 : size(missingGenes,1)
    Recon3DGeneInfo(r+1,'HGNC') = missingGenes(i,1);
    Recon3DGeneInfo(r+1,'EntrezGene') = missingGenes(i,2);
end

% Load translation of metabolite names from Metabolon to VMH from rBioNet
% (in MetaboAnnotator for the moment)
% this does not work for the moment as I used numeric ID and not metabolon
% name
%load('met_strc_rBioNet_new.mat');
%[VMH2IDmappingAll,VMH2IDmappingPresent,VMH2IDmappingMissing]=getIDfromMetStructure(metabolite_structure_rBioNet,'metabolon');
filenameMetabolon = 'metabolon_crossmatch_IT_withUpdatedInchiKey.xlsx';
[NUM,TXT,RAW]=xlsread(filenameMetabolon);
metabolon = RAW;
colMetabolon = find(contains(metabolon(1,:),'CHEMICAL_NAME'));
colVMH = find(contains(metabolon(1,:),'VMH'));
% append missing metabolon data
metabolon{end+1,colMetabolon} = 'thymine';
metabolon{end,colVMH} = 'thym';
metabolon{end+1,colMetabolon} = 'argininosuccinate';
metabolon{end,colVMH} = 'argsuc';
metabolon{end+1,colMetabolon} = '3-ureidoisobutyrate';
metabolon{end,colVMH} = '3uib';
metabolon{end+1,colMetabolon} = 'diacetylspermidine*';
metabolon{end,colVMH} = 'CE1059';
metabolon{end+1,colMetabolon} = 'carnosine';
metabolon{end,colVMH} = 'carn';
metabolon{end+1,colMetabolon} = 'homocarnosine';
metabolon{end,colVMH} = 'hmcarn';

% translate HNGC IDs in input file into EntrezGene IDs
colG = find(strcmp(lower(InputTable(1,:)),'gene_symbol'));

% translate metabolon name in input file into VMH IDs
colM = find(strcmp(lower(InputTable(1,:)),'biochemical'));
cnt = 0;

% find VMH col in dataset --DF
colVMHE = find(strcmp(lower(InputTable(1,:)),'vmhid')); 
clear geneMarkerListN
for i = 2 : size(InputTable,1)% first row is header
    gene = InputTable(i,colG);
    clear geneE
    for j = 1 : size(Recon3DGeneInfo,1)
        if strcmp(gene, Recon3DGeneInfo(j,colHGNC))
            geneE = (Recon3DGeneInfo{j,colE});
            geneH = (Recon3DGeneInfo{j,colHGNC});
        end
    end
    if exist('geneE','var')
        
        %met =  InputTable(i,colM);

        %for j = 1 : size(metabolon,1)
        %    if strcmp(met, metabolon(j,colMetabolon))
        %        M = (metabolon{j,colVMH});
        %    end
        %end
        
        clear M
        M = InputTable{i,colVMHE};
        
        % construct input into IEM modelling script
        %     geneMarkerList = {
        % %   %  gene  biomarker list
        %     '249.1' '3pg;cholp;glyc3p;ethamp'
        %     '1491.1' 'cyst_L' %?
        if  exist('M','var') && ~isempty(M) && ~strcmp(M, 'NaN') && ~strcmp(M, 'NA')
            if  ~isnan(M(1))
                geneEN=strcat(num2str(geneE), '.1');
                if ~exist('geneMarkerListN','var') || length(find(ismember(geneMarkerListN(:,1),geneEN)))==0
                    cnt = 1 + cnt';
                    geneMarkerListN{cnt,1} = geneEN;
                    geneMarkerListN{cnt,3} = geneH;

                end
                if size(geneMarkerListN,2)>1
                    
                    geneMarkerListN{ismember(geneMarkerListN(:,1),geneEN),2} = [geneMarkerListN{ismember(geneMarkerListN(:,1),geneEN),2},M,';'];
                else % add second col
                    geneMarkerListN{ismember(geneMarkerListN(:,1),geneEN),2} = [M,';'];
                end
            end
        end
    end
end
% now check that the metabolite is not in the list for the
% corresponding gene yet
clear geneMarkerListN2
for i = 1 : size(geneMarkerListN,1)
    clear M
    M = split(geneMarkerListN(i,2),';');
    M2 = unique(M);
    geneMarkerListN2{i,1} = geneMarkerListN{i,1};
    geneMarkerListN2{i,2} = M2{2};
    geneMarkerListN2{i,3} = geneMarkerListN{i,3};
 
    if length(M2)>2
        for j = 3 : length(M2)
            geneMarkerListN2{i,2} = strcat(geneMarkerListN2{i,2},';',M2{j});
        end
    end
end
geneMarkerList = geneMarkerListN2;
run_IEM_general;

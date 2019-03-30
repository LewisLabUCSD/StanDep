% % % % NOTE: You must have a working copy of COBRA Toolbox in MATLAB

% % % get data from Uhlen et al.
% rnaData = getUhlenData; % Source Uhlen et al., (ENSG)
% % % Alternatively, get cancer cell line data
% rnaData = getNCI60Data; % Source: CellMiner (Entrez)
% (OR)
% rnaData = load('data_RPKM_NCI60'); % Source: data used by Anne (Entrez)
% rnaData = rnaData.data;

% % load model with genes in ENSG format
load('ENSG_Recon22_Mapping','model'); % for tissue data
% % load model with genes in Entrez format
load('Entrez_Recon22_Mapping','model'); % for NCI60 data

% % get model data
modelData = getModelData(rnaData,model);

% % get enzyme data
% get specialist enzymes
spec = getSpecialistEnzymes(model);
% get promiscuous enzymes
prom = getPromEnzymes(model);
% compile into enzymatic data
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

% % generate cluster object
edgeX = [-3 -2 -1 0 1 2 2.5 3 4 5]; % choose this according the width of the distribution
k = 21; % chose this according to your needs
clustObjEnz = geneExprDist_hierarchy(enzymeData,[],edgeX,k,model,'euclidean','complete');

% % generate list of core reactions
rxnTisMat = models4mClusters1(clustObjEnz,enzymeData.Tissue,model,edgeX,[],[],false);

% % generate ubiquity matrix
ubiScore = getUbiquityScore(clustObjEnz,edgeX,model);

% % generate mcadre models
% write your code or function to do this
% tissueModel = mCADRE(model,ubiScore(:,i),rxnTisMat(:,i));
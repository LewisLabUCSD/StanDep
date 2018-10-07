function [rxnTisMat] = funcConstructModels(rnaData,model_file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% % load model with genes in ENSG format
load(model_file,'model');

%data_mat = 
modelData = getModelData(rnaData,model);

% % get enzyme data
% get specialist enzymes
spec = getSpecialistEnzymes(model);
% get promiscuous enzymes
prom = getPromEnzymes(model);
% compile into enzymatic data
enzymeData = comparePromiscuousSpecific(spec,prom,modelData,false);

% % generate cluster object
edgeX = [-3 -2 -1 0 1 2 2.5 3 4 5]; % choose this according the width of the distribution
k = 21; % chose this according to your needs
clustObjEnz = geneExprDist_hierarchy(enzymeData,[],edgeX,k,model,'euclidean','complete',false);

% % generate list of core reactions
rxnTisMat = models4mClusters1(clustObjEnz,enzymeData.Tissue,model,edgeX,'.',[],false,false);

% % generate ubiquity matrix
%ubiScore = getUbiquityScore(clustObjEnz,edgeX,model);


end


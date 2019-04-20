# StanDep
Codes for "StanDep: capturing transcriptomic variability improves context-specific metabolic models "
Chintan J. Joshi, Song-Min Schinn, Anne Richelle, Isaac Shamie, Eyleen J. O’Rourke, Nathan E. Lewis

# NOTES:
1. Do not add depreicated folder
2. You must have cobra toolbox on MATLAB

# Loading transcriptomics data (rnaData)
1. The transcriptomics data is required as a MATLAB structure containing following fields:
  a. gene: list of gene names
  b. value: a matrix whose each column contains the expression of genes across a tissue/context and each row contains expression of a gene across all tissues/contexts.
  c. genesymbols (optional): any alternative gene names that the user may want to keep track of.
  d. Tissue: names of conditions

2. exprData should contain only those entities that are in the model (genes/enzymes).
3. You must have access to the genome-scale model you want to use for integration.

# To get the active core reaction lists for the conditions follow the code below:
% % initialize COBRA toolbox

initCobraToolbox

% % extract expression data of the genes in the model

modelData = getModelData(rnaData,model);

% % calculate enzymes in the model

spec = getSpecialistEnzymes(model);
prom = getPromEnzymes(model);

% % calculate enzyme expression
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);

% % other inputs needed
edgeX = [-2 -1 0 1 2 2.5 3 4]; % bins
distMethod = 'euclidean'; % distance method
linkageMethod = 'complete'; % linkage metric for hierarchical clustering

% % calculate clusters of enzyme expression
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,k,distMethod,linkageMethod);

% % calculate active reaction lists as binary matrix
coreRxnMat = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],true,0,[1 1]); 

% % calculate Jaccard similarity between core reaction lists
calcJaccardSimilarity(coreRxnMat,enzymeData.Tissue,'matrix',true);

% % prepare model using mCADRE
ubiScore = getUbiquityScore(clustObj,edgeX,model); % calculate ubiquity score
% % % run mCADRE by using the above value for the input variable "ubiquityScore"
% % % Similarly other extraction methods can be tailored as well.

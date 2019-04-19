# StanDep
Codes for "StanDep: capturing transcriptomic variability improves context-specific metabolic models "
Chintan J. Joshi, Song-Min Schinn, Anne Richelle, Isaac Shamie, Eyleen J. O’Rourke, Nathan E. Lewis

# Loading transcriptomics data (rnaData)
1. The transcriptomics data is required as a MATLAB structure containing following fields:
  a. gene: list of gene names
  b. value: a matrix whose each column contains the expression of genes across a tissue/context and each row contains expression of a gene across all tissues/contexts.
  c. genesymbols (optional): any alternative gene names that the user may want to keep track of.
  d. Tissue: names of conditions

2. exprData should contain only those entities that are in the model (genes/enzymes).
3. You must have access to the genome-scale model you want to use for integration.

# To get the active core reaction lists for the conditions follow the code below:
initCobraToolbox
modelData = getModelData(rnaData,model);
spec = getSpecialistEnzymes(model);
prom = getPromEnzymes(model);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-2 -1 0 1 2 2.5 3 4];
distMethod = 'euclidean';
linkageMethod = 'complete';
clustObj = geneExprDist_hierarchy(exprData,[],edgeX,k,distMethod,linkageMethod)

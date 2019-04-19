function [rxnTisMat] = getLocalT1_case(modelData,model,localThs)

% USAGE:
% % [rxnTisMat] = getLocalT1_case(modelData,model,globalThs)

% INPUTS:
% % modelData:  gene expression data
                    % % value: a matrix of gene (rows) expression in conditions (cols)
                    % % gene: list of genes
                    % % Tissue: list of conditions
% % model:      a COBRA model where gene formats are same as modelData
% % localThs:  threshold value in percentile [0 100].

% OUTPUTS:
% % rxnTisMat:  a binary matrix indicating whether a reaction (row) is on a given condition (cols)

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)
% %                 ......implemented as Anne Richelle, please contact her

rawData = modelData.value;
linData = log10(reshape(rawData,numel(rawData),1));
linData(linData==-inf) = [];

if localThs <= 0
    error('Percentile threshold must be a positive value.');
else
    localThsVal1 = 10^(prctile(linData,localThs));
    maxGene = max(rawData,[],2);
    meanGene = mean(rawData,2);
    localGeneThs = maxGene;
    localGeneThs(meanGene>=localThsVal1) = meanGene(meanGene>=localThsVal1);
    % make selection of gene/enzyme-tissue pairs
    localSel1 = rawData>=localGeneThs;
    nLocalSel1 = sum(sum(rawData>=localGeneThs));
    
    spec = getSpecialistEnzymes(model);
    prom = getPromEnzymes(model);
    thrGeneData = modelData;
    thrGeneData.value = 5*log(1 + rawData./localGeneThs); % thresholding step
    thrEnzymeData = comparePromiscuousSpecific(spec,prom,thrGeneData); % convert gene data to enzyme data
    thrRxnData = convertEnzyme2RxnValues(thrEnzymeData,model); % convert enzyme data to reaction data
    rxnTisMat = thrRxnData.value > 5*log(2);
end
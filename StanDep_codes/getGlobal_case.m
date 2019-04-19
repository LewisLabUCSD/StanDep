function [rxnTisMat] = getGlobal_case(modelData,model,globalThs)

% USAGE:
% % [rxnTisMat] = getGlobal_case(modelData,model,globalThs)

% INPUTS:
% % modelData:  gene expression data
                    % % value: a matrix of gene (rows) expression in conditions (cols)
                    % % gene: list of genes
                    % % Tissue: list of conditions
% % model:      a COBRA model where gene formats are same as modelData
% % globalThs:  threshold value in percentile [0 100].

% OUTPUTS:
% % rxnTisMat:  a binary matrix indicating whether a reaction (row) is on a given condition (cols)

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)
% %                 ......implemented as Anne Richelle, please contact her

rawData = modelData.value;
linData = log10(reshape(rawData,numel(rawData),1));
linData(linData==-inf) = [];

if globalThs <= 0
    error('Percentile threshold must be a positive value.');
else
    n = size(rawData,1);
    globalThsVal = prctile(linData,globalThs);
    globalGeneThs = 10.^(repmat(globalThsVal,n,1));
    % make selection of gene-tissue pairs
    globalSel = rawData>=globalGeneThs;
    nGlobalSel = sum(sum(rawData>=globalGeneThs));
    
    spec = getSpecialistEnzymes(model);
    prom = getPromEnzymes(model);
    thrGeneData = modelData;
    thrGeneData.value = 5*log(1+rawData./globalGeneThs); % thresholding step
    thrEnzymeData = comparePromiscuousSpecific(spec,prom,thrGeneData); % convert gene data to enzyme data
    thrRxnData = convertEnzyme2RxnValues(thrEnzymeData,model); % convert enzyme data to reaction data
    rxnTisMat = thrRxnData.value > 5*log(2);
end
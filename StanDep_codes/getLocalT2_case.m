function [rxnTisMat] = getLocalT2_case(modelData,model,lowerThs,upperThs)

% USAGE:
% % [rxnTisMat] = getLocalT2_case(modelData,model,lowerThs,upperThs)

% INPUTS:
% % modelData:  gene expression data
                    % % value: a matrix of gene (rows) expression in conditions (cols)
                    % % gene: list of genes
                    % % Tissue: list of conditions
% % model:      a COBRA model where gene formats are same as modelData
% % lowerThs:   lower threshold value in percentile [0 100].
% % upperThs:   upper threshold value in percentile [0 100].

% OUTPUTS:
% % rxnTisMat:  a binary matrix indicating whether a reaction (row) is on a given condition (cols)

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)
% %                 ......implemented as Anne Richelle, please contact her

rawData = modelData.value;
linData = log10(reshape(rawData,numel(rawData),1));
linData(linData==-inf) = [];

if lowerThs > upperThs || lowerThs < 0 || upperThs < 0
    error('Thresholds are in percentiles, lower threshold should be less than upper threshold, and both thresholds should be positive.');
else
    % % calculate local thresholds T2
    localThsVal2 = 10^(prctile(linData,upperThs));
    localThsVal1 = 10^(prctile(linData,lowerThs));
    meanGene = mean(rawData,2);
    localGeneThs = meanGene;
    localGeneThs(meanGene>=localThsVal2) = localThsVal2;
    localGeneThs(meanGene<localThsVal1) = localThsVal1;
    % make selection of gene-tissue pairs
    localSel2 = rawData>=localGeneThs;
    nLocalSel2 = sum(sum(rawData>=localGeneThs));
    
    spec = getSpecialistEnzymes(model);
    prom = getPromEnzymes(model);
    thrGeneData = modelData;
    thrGeneData.value = 5*log(1 + rawData./localGeneThs); % thresholding step
    thrEnzymeData = comparePromiscuousSpecific(spec,prom,thrGeneData); % convert gene data to enzyme data
    thrRxnData = convertEnzyme2RxnValues(thrEnzymeData,model); % convert enzyme data to reaction data
    rxnTisMat = thrRxnData.value > 5*log(2);
end
function [rxnTisMat,globalSel,nGlobalSel] = getGlobal_case(modelData,model,globalThs)

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
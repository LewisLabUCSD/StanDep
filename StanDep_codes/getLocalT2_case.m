function [rxnTisMat,localSel2,nLocalSel2] = getLocalT2_case(modelData,model,lowerThs,upperThs)

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
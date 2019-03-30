function [rxnTisMat,localSel2,nLocalSel2] = getLocalT2(enzymeData,model,lowerThs,upperThs,tisid)

rawData = enzymeData.value;
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
    localGeneThs(meanGene<localThsVal1) = localThsVal1+1;
%     localGeneThs = sort(localGeneThs,'ascend');
    % make selection of gene-tissue pairs
    if tisid==0
        rxnTisMat = zeros(length(model.rxns),size(rawData,2));
        localSel2 = rawData>=localGeneThs;
        nLocalSel2 = sum(sum(rawData>=localGeneThs));
    else
        rxnTisMat = zeros(length(model.rxns),1);
        localSel2 = rawData(:,tisid)>=localGeneThs;
        nLocalSel2 = sum(sum(rawData(:,tisid)>=localGeneThs));
    end
    for i=1:size(localSel2,2)
        rxns_on = enzymeData.rxns(localSel2(:,i));
        rxns_on = unique([linearization_index(rxns_on(cellfun(@iscell,rxns_on)),'cols');rxns_on(~cellfun(@iscell,rxns_on))]);
        rxnTisMat(:,i) = ismember(model.rxns,rxns_on);
    end
end
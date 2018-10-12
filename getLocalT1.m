function [rxnTisMat,localSel1,nLocalSel1] = getLocalT1(enzymeData,model,localThs,tisid)

rawData = enzymeData.value;
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
    localGeneThs = sort(localGeneThs,'ascend');
    % make selection of gene/enzyme-tissue pairs
    if tisid==0
        rxnTisMat = zeros(length(model.rxns),size(rawData,2));
        localSel1 = rawData>=localGeneThs;
        nLocalSel1 = sum(sum(rawData>=localGeneThs));
    else
        rxnTisMat = zeros(length(model.rxns),1);
        localSel1 = rawData(:,tisid)>=localGeneThs;
        nLocalSel1 = sum(sum(rawData(:,tisid)>=localGeneThs));
    end
    for i=1:size(localSel1,2)
        rxns_on = enzymeData.rxns(localSel1(:,i));
        rxns_on = unique([linearization_index(rxns_on(cellfun(@iscell,rxns_on)),'cols');rxns_on(~cellfun(@iscell,rxns_on))]);
        rxnTisMat(:,i) = ismember(model.rxns,rxns_on);
    end
end
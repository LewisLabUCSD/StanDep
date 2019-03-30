function [rxnTisMat,globalSel,nGlobalSel] = getGlobal(enzymeData,model,globalThs,tisid)

rawData = enzymeData.value;
linData = log10(reshape(rawData,numel(rawData),1));
linData(linData==-inf) = [];

if globalThs <= 0
    error('Percentile threshold must be a positive value.');
else
    n = size(rawData,1);
    globalThsVal = prctile(linData,globalThs);
    globalGeneThs = 10.^(repmat(globalThsVal,n,1));
    % make selection of gene-tissue pairs
    if tisid==0
        rxnTisMat = zeros(length(model.rxns),size(rawData,2));
        globalSel = rawData>=globalGeneThs;
        nGlobalSel = sum(sum(rawData>=globalGeneThs));
    else
        rxnTisMat = zeros(length(model.rxns),1);
        globalSel = rawData(:,tisid)>=globalGeneThs;
        nGlobalSel = sum(sum(rawData(:,tisid)>=globalGeneThs));
    end
    for i=1:size(globalSel,2)
        rxns_on = enzymeData.rxns(globalSel(:,i));
        rxns_on = unique([linearization_index(rxns_on(cellfun(@iscell,rxns_on)),'cols');rxns_on(~cellfun(@iscell,rxns_on))]);
        rxnTisMat(:,i) = ismember(model.rxns,rxns_on);
    end
end
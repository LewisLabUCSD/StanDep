function [modelData] = getModelData(expressionData,model)

% if gene ids are in double, convert them to string
if isnumeric(expressionData.gene)
    expressionData.gene = cellstr(char(num2str(expressionData.gene)));
    expressionData.gene = strrep(expressionData.gene,' ','');
end
geneNames = expressionData.gene; 

% construct model data structure
if isfield(expressionData,'symbol')
    modelIdx = ismember(expressionData.gene,model.genes) | ismember(expressionData.symbol,model.genes);
else
    modelIdx = ismember(expressionData.gene,model.genes);
end
modelData.gene = geneNames(modelIdx);
if isfield(expressionData,'meanValue')
    modelData.meanValue = expressionData.meanValue(modelIdx);
end
if isfield(expressionData,'unit')
    modelData.unit = expressionData.unit;
end
if isfield(expressionData,'Tissue')
    modelData.Tissue = expressionData.Tissue;
end
if isfield(expressionData,'genesymbol')
    modelData.genesymbol = expressionData.genesymbol(modelIdx);
end
for j=1:size(expressionData.valuebyTissue,2)
    modelData.value(:,j) = expressionData.valuebyTissue(modelIdx,j);
end

% sore missing gene IDs
modelData.ID_geneMissing = model.genes(~ismember(model.genes,expressionData.gene));
modelData.ID_genePresent = model.genes(ismember(model.genes,expressionData.gene));

% resolve non-uniques
geneids = unique(modelData.gene);
data = zeros(length(geneids),length(modelData.Tissue));
for i=1:length(geneids)
    ix = ismember(modelData.gene,geneids{i});
    if sum(ix)>1
        data(i,:) = mean(modelData.value(ix,:),1);
    else
        data(i,:) = modelData.value(ix,:);
    end
end
modelData.value = data;
modelData.gene = geneids;
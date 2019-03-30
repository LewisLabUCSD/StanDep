function [enzymeData] = comparePromiscuousSpecific(spec,prom,modelData)

% comparing distributions of specialist & promiscuous enzymes
specSubunits = regexp(spec.enzymes,' & ','split');
promSubunits = regexp(prom.enzymes,' & ','split');

specExprMatrix = zeros(size(specSubunits,1),size(modelData.value,2));
for i=1:size(modelData.value,2)
    for j=1:size(specSubunits,1)
        if sum(ismember(modelData.gene,specSubunits{j}))==length(specSubunits{j})
            specExprMatrix(j,i) = min(modelData.value(ismember(modelData.gene,specSubunits{j}),i));
        end
    end
end
promExprMatrix = zeros(size(promSubunits,1),size(modelData.value,2));
for i=1:size(modelData.value,2)
    for j=1:size(promSubunits)
        if sum(ismember(modelData.gene,promSubunits{j}))==length(promSubunits{j})
            promExprMatrix(j,i) = min(modelData.value(ismember(modelData.gene,promSubunits{j}),i));
        end
    end
end
enzAbund = [specExprMatrix; promExprMatrix];
specExprVec = reshape(specExprMatrix,numel(specExprMatrix),1);
promExprVec = reshape(promExprMatrix,numel(promExprMatrix),1);

specExprVec(specExprVec==-1) = [];
promExprVec(promExprVec==-1) = [];

figure;
hold on
histogram(log10(specExprVec),'BinEdges',[-8:0.25:5],'DisplayStyle','stairs');
histogram(log10(promExprVec),'BinEdges',[-8:0.25:5],'DisplayStyle','stairs');
hold off
legend('Specialist Enzymes','Promiscuous Enzymes','Location','NorthWest');

enzymeData.enzyme = [spec.enzymes;prom.enzymes];
enzymeData.value = enzAbund;
enzymeData.rxns = [spec.rxns;prom.rxns];
enzymeData.Tissue = modelData.Tissue;
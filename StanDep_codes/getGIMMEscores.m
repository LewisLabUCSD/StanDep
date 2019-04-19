function gimmeScores = getGIMMEscores(clustObj,edgeX,model)

% USAGE:
% % gimmeScores = getGIMMEscores(clustObj,edgeX,model)
% % code needed to calculate inputs for GIMME

% INPUTS:
% % clustObj:   cluster object calculated in geneExprDist_hierarchy
% % edgeX:      bins used in clustObj
% % model:      a COBRA model to be used

% OUTPUTS:
% % gimmeScore: a matrix describing scores for each reaction

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

objDist = zeros(size(clustObj.Data));
uci = 1:1:size(clustObj.C,1);
cidx = clustObj.cindex;
[~,thrVal] = clusterVariability1(clustObj,edgeX,false,0,[1 1]);
md = max(max(clustObj.Data - thrVal(clustObj.cindex)'));
for i=1:length(uci)
    ic = find(cidx==uci(i));
    for j=1:size(clustObj.Data,2)
        objDist(ic,j) = clustObj.Data(ic,j) - thrVal(i);
    end
    objDist(ic,:) = objDist(ic,:)*md/max(max(objDist(ic,:)));
end
gimmeScores = objDist;
gimmeScores = transformToModel(gimmeScores,clustObj,model);

function B = transformToModel(A,clustObj,model)

[A,rxns] = openWeightMatrix(clustObj,A);
B = repmat(-1,length(model.rxns),size(clustObj.Data,2));
for i=1:length(model.rxns)
    if sum(ismember(rxns,model.rxns{i}))~=0
        B(i,:) = max(A(ismember(rxns,model.rxns{i}),:),[],1);
    end
end

function [openUbi,rxns] = openWeightMatrix(clustObj,uMat)

k = 0;
for i=1:size(uMat,1)
    if ischar(clustObj.objectMaps{i})
        k = k+1;
        openUbi(k,:) = uMat(i,:);
        rxns{k,1} = clustObj.objectMaps{i};
    else
        n = length(clustObj.objectMaps{i});
        for j=1:n
            k = k+1;
            openUbi(k,:) = uMat(i,:);
            rxns{k,1} = clustObj.objectMaps{i}{j,1};
        end
    end        
end
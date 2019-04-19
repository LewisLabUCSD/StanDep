function [H,M] = getMBAsets(clustObj,edgeX,model,tol)

% USAGE:
% % [H,M] = getMBAsets(clustObj,edgeX,model,tol)
% % code needed to calculate inputs for MBA

% INPUTS:
% % clustObj:   cluster object calculated in geneExprDist_hierarchy
% % edgeX:      bins used in clustObj
% % model:      a COBRA model to be used
% % tol:        a tolerance within which medium confidence reactions are placed

% OUTPUTS:
% % H:          a matrix describing whether a reaction is high confidence
% % M:          a matrix describing whether a reaction is medium confidence

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

objDist_lb = zeros(size(clustObj.Data));
objDist_ub = objDist_lb;
H = objDist_lb;
M = objDist_lb;
uci = 1:1:size(clustObj.C,1);
cidx = clustObj.cindex;
[~,thrVal] = clusterVariability1(clustObj,edgeX,false,0,[1 1]);

for j=1:size(clustObj.Data,2)
    for i=1:length(uci)
        ic = find(cidx==uci(i));
        objDist_ub(ic,j) = clustObj.Data(ic,j) - thrVal(i)*(1+tol);
        objDist_lb(ic,j) = clustObj.Data(ic,j) - thrVal(i)*(1-tol);
    end
end
H(objDist_ub > 0) = 1;
M(objDist_lb > 0 & objDist_ub < 0) = 1;

H = transformToModel(H,clustObj,model);
M = transformToModel(M,clustObj,model);
M(H==1 & M==1) = 0;

function B = transformToModel(A,clustObj,model)

[A,rxns] = openUbiquityMatrix(clustObj,A);
B = zeros(length(model.rxns),size(clustObj.Data,2));
for i=1:length(model.rxns)
    if sum(ismember(rxns,model.rxns{i}))~=0
        B(i,:) = max(A(ismember(rxns,model.rxns{i}),:),[],1);
    end
end

function [D,rxns] = openUbiquityMatrix(clustObj,A)

k = 0;
for i=1:size(A,1)
    if ischar(clustObj.objectMaps{i})
        k = k+1;
        D(k,:) = A(i,:);
        rxns{k,1} = clustObj.objectMaps{i};
    else
        n = length(clustObj.objectMaps{i});
        for j=1:n
            k = k+1;
            D(k,:) = A(i,:);
            rxns{k,1} = clustObj.objectMaps{i}{j,1};
        end
    end        
end

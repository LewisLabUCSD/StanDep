function [J,comparisonMade] = evaluateClusterNumber_NCI60(nc)

% evauates the influence of number of clusters
% mean Jaccard similarity across cell lines and cluster numbers = 0.9090


% get enzyme data for NCI60
load('E:\chintan_NCI60_consistentGEMmodels\NCI60Clusters\rxn_clusters\nci60_enzyme_fastcore_models','clustObj','edgeX','enzymeData','model');
k = 0;
nc = nc:1:max(clustObj.cindex);
comparisonMade = zeros(length(nc)*(length(nc)-1)/2,2);
for i=1:length(nc)
    fprintf('Running for %d clusters...\n',nc(i));
    newClustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,nc(i),model,'euclidean','complete');
    rxnTisMat{i} = models4mClusters1(newClustObj,enzymeData.Tissue,model,edgeX,[],[],false,0);
    Nr(:,i) = sum(rxnTisMat{i},1);
    close all
    if i>1
        for j=1:i-1
            k = k+1;
            J(:,k) = compareTissuesForTwoSets(rxnTisMat{i},rxnTisMat{j});
            comparisonMade(k,1) = nc(i);
            comparisonMade(k,2) = nc(j);
        end
    end
end

function Jacc = compareTissuesForTwoSets(set1,set2)

Jacc = zeros(size(set1,2),1);
for i=1:size(set1,2)
    Jacc(i) = sum(set1(:,i) & set2(:,i))/sum(set1(:,i) | set2(:,i));
end
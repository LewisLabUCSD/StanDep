function [geneTis,cutOff,thr] = models4mClusters1(clustObj,tisNames,model,xedge,folderName,cutOff,figFlag)

if isempty(cutOff)
    cutOff = clusterVariability1(clustObj,xedge,figFlag);
%     data = reshape(clustObj.Data,numel(clustObj.Data),1); data(data==-inf) = 0;
%     thr = prctile(data,cutOff);
end
nTis = size(clustObj.Data,2);
% % zscore the entire data
% nGenes = size(clustObj.Data,1);
% data = reshape(clustObj.Data,numel(clustObj.Data),1);
% dmu = nanmean(data); dsig = nanstd(data);
% data = (data-repmat(dmu,length(data),1))./repmat(dsig,length(data),1);
% clutObj.Data = reshape(data,nGenes,nTis);

nClust = size(clustObj.C,1);
tisNames = strrep(tisNames,'_',' ');
geneTis = false(length(clustObj.objects),nTis);
for i=1:nClust
    gid = clustObj.cindex==i;
% % % % %     when threholds are applied for each tissue cluster pair
% % % %     gl = clustObj.objects(gid);
% % % %     clustMat = false(sum(gid),nTis);
% % % %     for j=1:nTis
% % % %         clustData = clustObj.Data(gid,i);
% % % %         thr = prctile(clustData,roundn(cutOff(i),-2));
% % % %         ind = clustData >= thr;
% % % %         g = gl(ind);
% % % %         geneTis(:,j) = geneTis(:,j) | ismember(clustObj.objects,g);
% % % %         clustMat(:,j) = ismember(clustObj.objects(gid),g);
% % % %     end
    % when thresholds are applied for each cluster globally
    clustMat = false(sum(gid),nTis);
    clustData = clustObj.Data(gid,:);
%     clustData = clustData./max(clustData,[],2);
    clustData = reshape(clustData,numel(clustData),1);
    colgene = repmat(clustObj.objects(gid),nTis,1);
    coltis = repmat(1:1:nTis,sum(gid),1);
    coltis = reshape(coltis,numel(coltis),1);
    remi = clustData==-inf;
    coltis(remi) = []; colgene(remi) = [];
    clustData(remi) = [];
    thr(i,1) = prctile(clustData,roundn(cutOff(i),-2));
    ind = clustData >= thr(i);
    fprintf('fraction selected = %0.4f\n',sum(ind)/(sum(gid)*nTis));
    coltis = coltis(ind);
    colgene = colgene(ind);
    for j=1:nTis
        if sum(coltis==j)~=0
            geneTis(:,j) = geneTis(:,j) | ismember(clustObj.objects,unique(colgene(coltis==j)));
            clustMat(:,j) = ismember(clustObj.objects(gid),unique(colgene(coltis==j)));
        end
    end
%     % when thresholds are applied as local T1 for each cluster
%     gl = clustObj.objects(gid);
%     clustMat = false(sum(gid),nTis);
%     clustData = clustObj.Data(gid,:);
%     linClust = reshape(clustData,numel(clustData),1);
%     linClust(linClust==-inf) = [];
%     thr = prctile(linClust,roundn(cutOff(i),-2));
%     for j=1:size(clustData,1)
%         if mean(clustData(j,:)) >= thr
%             clustMat(j,:) = true;
%             geneTis(strcmp(clustObj.objects,gl(j)),:) = true;
%         else
%             clustMat(j,:) = clustData(j,:) > thr;
%             geneTis(strcmp(clustObj.objects,gl(j)),:) = clustData(j,:) > thr;
%         end
%         if sum(clustData(j,:)==-inf)~=0
%             clustMat(j,clustData(j,:)==-inf) = false;
%             geneTis(strcmp(clustObj.objects,gl(j)),clustData(j,:)==-inf) = false;
%         end
%     end
    if ~isempty(folderName)
        Jc = getJaccardSimMatrix(clustMat);
        writeMyMatrix(Jc,strcat(folderName,'/C',num2str(i),'.txt'));
        if i==1
            Z = linkage(Jc,'complete');
            Y = pdist(Jc,'euclidean');
            leafOrderg = optimalleaforder(Z,Y);
            fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
                'complete','euclidean',cophenet(Z,Y));
            if figFlag
                figure;
                [Hg,Tg,outpermg] = dendrogram(Z,0,'Orientation','left','labels',tisNames,'Reorder',leafOrderg);
                title('Tissues clustered based on metabolic similarity in Cluster 1');
            end
        end
    end
end

J = getJaccardSimMatrix(geneTis);
writeMyMatrix(J,strcat(folderName,'/Call.txt'));

if figFlag
    Z = linkage(J,'complete');
    Y = pdist(J,'euclidean');
    leafOrderg = optimalleaforder(Z,Y);
    fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
        'complete','euclidean',cophenet(Z,Y));
    figure;
    [Hg,Tg,outpermg] = dendrogram(Z,0,'Orientation','left','labels',tisNames,'Reorder',leafOrderg);
    title('Tissues clustered based on metabolic similarity');
    J = J(outpermg,outpermg);
    tisNames = tisNames(outpermg);

    figure;
    imagesc(J);
    ax = gca; ax.TickDir = 'out'; ax.DataAspectRatioMode = 'manual';
    xticks(1:1:nTis); yticks(1:1:nTis);
    xticklabels(tisNames);
    yticklabels(tisNames);
    xtickangle(45);
    if ~isfield(clustObj,'objectMaps')
        title('Jaccard Similarity for genes in tissues');
    else
        title('Jaccard Similarity for enzymes in tissues');
    end
    cbi = colorbar;
    cbi.TickDirection = 'out';
    colormap jet
end

if ~isempty(model)
    if ~isfield(clustObj,'objectMaps')
        model.genes = unique(model.genes);
        activeGenes = cell(1,nTis);
        for i=1:nTis
            activeGenes{i} = find(ismember(model.genes,clustObj.objects(geneTis(:,i))));
        end
        geneTis = activeGenes;
    else
        activeRxns = false(length(model.rxns),nTis);
        for i=1:nTis
            mylist = clustObj.objectMaps(geneTis(:,i));
            l1 = mylist(~cellfun(@iscell,mylist));
            l2 = linearization_index(mylist(cellfun(@iscell,mylist)),'cols');
            mylist = unique([l1; l2]);
            activeRxns(:,i) = ismember(model.rxns,mylist);
        end
        geneTis = activeRxns;
    end
end


function J = getJaccardSimMatrix(A)

n = size(A,2);
J = zeros(n,n);
for i=1:n
    for j=i:n
        J(j,i) = sum(A(:,i) & A(:,j))/sum(A(:,i) | A(:,j));
        J(i,j) = sum(A(:,i) & A(:,j))/sum(A(:,i) | A(:,j));
    end
end

function writeMyMatrix(A,fileName)

fid = fopen(fileName,'wt');

for ii = 1:size(A,1)
    fprintf(fid,'%0.2f\t',A(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);
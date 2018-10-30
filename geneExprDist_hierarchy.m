function [clustObj,Z,Hg,Tg,outpermg,tisMat,rowNames] = geneExprDist_hierarchy(exprData,removeObjects,edgeX,k,model,distMethod,linkageMethod,varargin)

%% See if figure should be plotted or not
numvarargs = length(varargin);
if numvarargs > 1
    error('TooManyInputs', ...
        'requires at most 1 optional inputs');
end
% set defaults for optional inputs
optargs = {true};
optargs(1:numvarargs) = varargin;
% Place optional args as variable names
[to_plot] = optargs{:};

%% 
if isfield(exprData,'value')
    exprData.valuebyTissue = exprData.value;
    clear exprData.value
end

if ~isempty(removeObjects)
    remIndex = find(ismember(exprData.gene,removeObjects));
    exprData.gene(remIndex) = [];
    exprData.valuebyTissue(remIndex,:) = [];
end
if  isfield(exprData,'gene')
    exprData.gene(sum(exprData.valuebyTissue,2)==0) = [];
elseif isfield(exprData,'enzyme')
    enz=1;
    exprData.enzyme(sum(exprData.valuebyTissue,2)==0) = [];
    exprData.rxns(sum(exprData.valuebyTissue,2)==0) = [];
end
exprData.valuebyTissue(sum(exprData.valuebyTissue,2)==0,:) = [];
nGenes = size(exprData.valuebyTissue,1);
nTis = size(exprData.valuebyTissue,2);

% v = exprData.valuebyTissue./max(exprData.valuebyTissue,[],1);
% v = log10(v);

v = log10(exprData.valuebyTissue);

% v = exprData.valuebyTissue*max(max(exprData.valuebyTissue,[],1))./max(exprData.valuebyTissue,[],1);
% v = log10(v); 

% v = plotDoubleNormalizedData(exprData,false);

% figure out number of number of columns and rows in figure
nrow = 2;
ncol = 2;
while nrow*ncol <= k
    ncol = ncol + 1;
    if nrow*ncol <= k
        nrow = nrow + 1;
    end
end

Yhist = zeros(nGenes,length(edgeX)-1);
for i=1:nGenes
    temp = v(i,:);
    temp(temp==-inf) = [];
    [Yhist(i,:)] = histcounts(temp,edgeX);
end
Yhist = Yhist/nTis; % normalize by total number of tissues
% gTis = sum(Yhist,2);
% Yhist = Yhist./max(Yhist,2); % normalize each gene distribution by max

% hierarchical clustering
Y = pdist(Yhist,distMethod);
% y = squareform(Y);
Z = linkage(Y,linkageMethod);
% Y = pdist(Yhist,distMethod);
leafOrderg = optimalleaforder(Z,Y);

%fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
%    linkageMethod,distMethod,cophenet(Z,Y));

% idx = cluster(Z,'maxclust',k,'criterion','distance');
% cObj = clustergram(Yhist);
if to_plot
    figure;
else
    figure('visible','off');
end

[Hg,Tg,outpermg] = dendrogram(Z,k,'Orientation','left','Reorder',leafOrderg);
idxg = Tg;

Cg = zeros(k,length(edgeX)-1);
for i=1:k
    Cg(i,:) = mean(Yhist(idxg==i,:),1);
end
clustObj.Distribution = Yhist;
numObjs = roundn(sum(Yhist,2).*nTis,0);
% numObjs = roundn(sum(Yhist,2).*gTis,0);
clustObj.numObjs = numObjs;
if  isfield(exprData,'gene')
    clustObj.objects = exprData.gene;
elseif isfield(exprData,'enzyme')
    clustObj.objects = exprData.enzyme;
    clustObj.objectMaps = exprData.rxns;
end
if isfield(exprData,'genesymbol')
    clustObj.altObjects = exprData.genesymbol;
end
clustObj.C = Cg;
clustObj.cindex = idxg;
clustObj.Data = v;
[idxg,ig] = sort(idxg,'ascend');
Yhist = Yhist(ig,:); 
% v = v(ig,:);
% clustObj.objects = exprData.gene(ig);

% for i=1:size(C,1)
%     ia = (idx == i);
%     if i==1
%         sepindex(i,1) = sum(ia);
%     else
%         sepindex(i,1) = sum(ia) + sepindex(i-1,1);
%     end
% end

xtiklab = cell(length(edgeX)-1,1);
for i=1:length(edgeX)-1
    xtiklab{i} = ['[' num2str(edgeX(i)) ',' num2str(edgeX(i+1)) ']'];
end
% figure;
% imagesc(Yhist*nTis); hold on
% ax = gca; ax.TickDir = 'out';
% % imagesc(Yhist.*gTis); hold on
% % plot([0.5 length(edgeX)-.5],[sepindex sepindex],'Color',[0.74 0.26 0.1],'LineWidth',2);
% xticks(1:1:(length(edgeX)-1));
% xticklabels(xtiklab);
% xlabel('log_{10}[expr.]');
% ylabel('gene #');
% hold off

% % colormap(cubehelix([],0.1,0.4,2,1,[0.1,1],[0,0.9]));
% colormap(cubehelix([],0.5,-1.5,1,1,[0.1,1],[0 0.9]));
% cbv = colorbar('v');
% cbv.Ticks = 4:4:nTis;
% cbv.Label.String = '# of tissues';
% cbv.Label.FontSize = 14;
% cbv.TickDirection = 'out';

% if ~isempty(model)
%     figure;
%     trRxns = findTransRxns(model);
%     tr.genes = model.genes(sum(model.rxnGeneMat(ismember(model.rxns,trRxns),:),1)~=0);
%     ntr.genes = model.genes(sum(model.rxnGeneMat(~ismember(model.rxns,trRxns),:),1)~=0);
%     % all model genes
%     getModelClusters(clustObj,model,[0.4 0 0.4]); hold on
%     % non-transporter genes
%     getModelClusters(clustObj,ntr,[0.4 0 0]);
%     % transporter genes
%     getModelClusters(clustObj,tr,[0 0.4 0]); hold off
%     legend('all genes','metabolic genes','transporter genes');
% end

% figure;
% silhouette(Yhist,idx,distMethod);

if to_plot
    figure;
else
    figure('visible','off');
end
cbg = cubehelix(nTis+1,0.6,-0.9,1.9,1.1,[0.2 0.8],[0.4 0.8]);
for i=1:k
    subplot(nrow,ncol,i);
    gid = idxg==i;
%     h = plot(1:1:length(xtiklab),Yhist(gid,:).*gTis(gid)); hold on
    h = plot(1:1:length(xtiklab),Yhist(gid,:)'.*nTis); hold on
    ax = gca; ax.TickDir = 'out';
    set(h,{'Color'},num2cell(cbg(nTis+1-numObjs(gid),:),2));
    if i==k
        colormap(cubehelix([],0.6,-0.9,1.9,1.1,[0.2 0.8],[0.4 0.8]));
        cbi = colorbar('Position',[0.925 0.18 0.02 0.7]);
        cbi.Ticks = linspace(0,1,5);
        cbi.TickLabels = num2cell(nTis:-nTis/4:0);
        cbi.Label.String = 'Number of tissues gene is expressed in';
        cbi.Label.FontSize = 14;
        cbi.TickDirection = 'out';
    end
    plot(1:1:length(xtiklab),Cg(i,:)*nTis,'Color','k','LineWidth',2); hold off
    xlim([1 (length(edgeX)-1)]);
    xticks(1:1:(length(edgeX)-1));
    xticklabels(xtiklab);
    xlabel('Bins of log_{10}[expr.]');
    xtickangle(45);
    clustObj.numObjInClust(i,1) = sum(gid);
    if isfield(exprData,'gene')
        title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' genes']);
    elseif isfield(exprData,'enzyme')
        title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' enzymes']);
    end
    ylim([0 nTis]);
    yticks(0:0.25*nTis:nTis);
    yticklabels(0:0.25*nTis:nTis);
    ylabel('# of tissues');
end

if to_plot
    figure;
else
    figure('visible','off');
end
tisMat = []; % initialize tissueMatrix
rowNames = [];
r = mod(nTis,4);
% Q = fix(nTis,4);
for i=1:k
    gid = clustObj.cindex==i;
    ngenes = sum(gid);
    geneTisMat = zeros(nTis,length(xtiklab));
    for j=1:nTis
        temp = v(gid,:);
        temp = temp(:,j);
        temp(temp==-inf) = [];
        geneTisMat(j,:) = histcounts(temp,edgeX)';
    end;
    tisMat = [tisMat;geneTisMat/ngenes]; % join tissueMatrix
    myStrings = strcat(exprData.Tissue,'_',num2str(i));
    rowNames = [rowNames;myStrings];
    geneTisMat = 1-geneTisMat/ngenes;
    subplot(nrow,ncol,i);
    imagesc(geneTisMat);
    caxis([0 1]);
    ax = gca; ax.TickDir = 'out';
    ylim([0 nTis]);
    yticks(0:0.25*(nTis-r):(nTis-r));
    yticklabels(0:0.25*(nTis-r):(nTis-r));
    ylabel('tissue #');
    xlim([0.5 (length(edgeX)-0.5)]);
    xticks(1:1:(length(edgeX)-1));
    xticklabels(xtiklab);
    xlabel('Bins of log_{10}[expr.]');
    xtickangle(45);
    colormap(cubehelix([],0.1,0.4,2,1,[0.1 1],[0 1]));
    if isfield(exprData,'gene')
        title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' genes']);
    elseif isfield(exprData,'enzyme')
        title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' enzymes']);
    end
end
colormap(cubehelix([],0.1,0.4,2,1,[0.1 1],[0 1]));
cbi = colorbar('Position',[0.925 0.18 0.02 0.7]);
cbi.Ticks = linspace(0,1,6);
cbi.TickLabels = num2cell(1:-1/5:0);
cbi.Label.String = 'fraction of genes in the cluster';
cbi.Label.FontSize = 14;
cbi.TickDirection = 'out';

% r = mod(k,8);
% Q = fix(k/8);
% if r==0
%     outloops = Q;
% else
%     outloops = Q+1;
% end
% 
% for f=1:outloops
%     figure;
%     if f~=outloops
%         for i = (f-1)*8+1:f*8
%             gid = clustObj.cindex==i;
%             ngenes = sum(gid);
%             gtMat = zeros(nTis,length(xtiklab));
%             for j=1:nTis
%                 temp = v(gid,:);
%                 temp = temp(:,j);
%                 temp(temp==-inf) = [];
%                 gtMat(j,:) = histcounts(temp,edgeX)';
%             end
%             gtMat = gtMat/ngenes;
%             Z = linkage(gtMat,linkageMethod);
%             Y = pdist(gtMat,distMethod);
%             leafOrderg = optimalleaforder(Z,Y);
%             fprintf('Cluster # %d: Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
%                 i,linkageMethod,distMethod,cophenet(Z,Y));
%             subplot(2,4,i-(f-1)*8);
%             dendrogram(Z,0,'Orientation','left','Labels',strrep(exprData.Tissue,'_','-'),'Reorder',leafOrderg);
%             ax = gca; ax.TickDir = 'out';
%             title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' genes']);
%         end
%     else
%         for i = (f-1)*8+1:k
%             gid = clustObj.cindex==i;
%             ngenes = sum(gid);
%             gtMat = zeros(nTis,length(xtiklab));
%             for j=1:nTis
%                 temp = v(gid,:);
%                 temp = temp(:,j);
%                 temp(temp==-inf) = [];
%                 gtMat(j,:) = histcounts(temp,edgeX)';
%             end
%             gtMat = gtMat/ngenes;
%             Z = linkage(gtMat,linkageMethod);
%             Y = pdist(gtMat,distMethod);
%             leafOrderg = optimalleaforder(Z,Y);
%             fprintf('Cluster # %d: Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
%                 i,linkageMethod,distMethod,cophenet(Z,Y));
%             subplot(2,4,i-(f-1)*8);
%             dendrogram(Z,0,'Orientation','left','Labels',strrep(exprData.Tissue,'_','-'),'Reorder',leafOrderg);
%             ax = gca; ax.TickDir = 'out';
%             title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' genes']);
%         end
%     end
% end

% figure;
% xedge = -2:0.5:5;
% xtiklab_ = cell(length(xedge)-1,1);
% for i=1:length(xedge)-1
%     xtiklab_{i} = ['[' num2str(xedge(i)) ',' num2str(xedge(i+1)) ']'];
% end
% geneTisMat = zeros(k,(length(xedge)-1));
% for i=1:k
%     gid = clustObj.cindex==i;
%     temp = v(gid,:);
%     temp(temp==-inf) = [];
%     geneTisMat(i,:) = histcounts(temp,xedge,'normalization','pdf')';
%     subplot(nrow,ncol,i);
%     bar(1:1:length(xtiklab_),geneTisMat(i,:),'BarWidth',1,'EdgeColor','none');
%     ax = gca; ax.TickDir = 'out';
%     ylim([0 1]);
%     yticks(0:0.25:1);
% %     yticklabels(0:0.25*nTis:nTis);
%     ylabel('# of tissue-genes');
%     xlim([0.5 (length(xedge)-0.5)]);
%     xticks(1:1:length(xedge));
%     xticklabels(xtiklab_);
%     xlabel('Bins of log_{10}[expr.]');
%     xtickangle(45);
%     set(gca,'fontsize',8);
% %     colormap(cubehelix([],0.1,0.4,2,1,[0.1 1],[0 1]));
%     if isfield(exprData,'gene')
%         title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' genes']);
%     elseif isfield(exprData,'enzyme')
%         title(['Cluster # ' num2str(i) '; ' num2str(sum(gid)) ' enzymes']);
%     end
% end

% % cluster tissueMatrix hierarchically
% % hierarchical clustering
% % Y = pdist(Yhist,distMethod);
% % y = squareform(Y);
% Zt = linkage(tisMat,linkageMethod);
% Yt = pdist(tisMat,distMethod);
% leafOrdert = optimalleaforder(Zt,Yt);
% fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
%     linkageMethod,distMethod,cophenet(Zt,Yt));
% % idx = cluster(Z,'maxclust',k,'criterion','distance');
% % cObj = clustergram(Yhist);
% figure;
% [Ht,Tt,outpermt] = dendrogram(Zt,'ColorThreshold',0.1,'Orientation','left','Reorder',leafOrdert);
% ax = gca; ax.TickDir = 'out';
% idxt = Tt;
% 
% % Ct = zeros(37,length(edgeX)-1);
% % for i=1:k
% %     Ct(i,:) = mean(tisMat(idxt==i,:),1);
% % end
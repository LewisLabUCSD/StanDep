function [cutOff,thrval,m,s,term1,term2,stdExpr,meanExpr] = clusterVariability1(clustObj,xedge,figFlag)

% % zscore the entire data
% nTis = size(clustObj.Data,2);
% nGenes = size(clustObj.Data,1);
% data = reshape(clustObj.Data,numel(clustObj.Data),1);
% dmu = nanmean(data); dsig = nanstd(data);
% data = (data-repmat(dmu,length(data),1))./repmat(dsig,length(data),1);
% clutObj.Data = reshape(data,nGenes,nTis);

cidx = clustObj.cindex;
uci = 1:1:size(clustObj.C,1);

meanExpr = zeros(size(cidx,1),1);
stdExpr = zeros(size(cidx,1),1);
varExpr = zeros(size(cidx,1),1);
noiseExpr = zeros(size(cidx,1),1);
for i=1:length(cidx)
    tempo = clustObj.Data(i,:);
    tempo(tempo==-inf) = [];
    meanExpr(i) = mean(tempo); % gene-specific mean
    stdExpr(i) = std(tempo); % gene-specific variability
    varExpr(i) = var(tempo);
    noiseExpr(i) = stdExpr(i)/meanExpr(i);
end
muData = clustObj.Data; muData(muData==-inf) = [];
muData = reshape(muData,numel(muData),1);
fprintf('Top 25th percentile for the data = %0.4f\n',prctile(muData,75));
muData = mean(muData);
fprintf('Mean of Data = %0.4f\n',muData);
sigData = sqrt(mean(varExpr));
fprintf('Std. Dev. of Data = %0.4f\n',sigData);
s = zeros(length(uci),1); m = s; mu_sm = s; noi = s;
if figFlag
    figure;
    subplot(2,2,4)
end
for i=1:length(uci)
    % variability in the mean of the genes in that cluster
    s(i) = sqrt(mean(varExpr(cidx==uci(i))));
    m(i) = mean(meanExpr(cidx==uci(i)));
    mu_sm(i) = mean(stdExpr(cidx==uci(i))./meanExpr(cidx==uci(i)));
    noi(i) = mean(noiseExpr(cidx==uci(i)));
end
if figFlag
    histogram(stdExpr);
    xlabel('\sigma_{[FPKM]},log_e');
    ylabel('# of genes');
    
    subplot(2,2,1)
    hold on
    for i=1:length(uci)
        plot(i,stdExpr(cidx==uci(i)),...
            'Marker','o','MarkerSize',4,'markerfacecolor','k','markeredgecolor','none');
    end
    plot(uci,s,...
            'Marker','o','MarkerSize',8,'markerfacecolor','r','markeredgecolor','none','linestyle','none');
    plot(uci,repmat(sigData,size(uci)),'--');
    hold off
    ylabel('\mu_{\sigma} of cluster');
    xlabel('cluster #');
    xlim([0 length(uci)+1]);
    
    subplot(2,2,2)
    hold on
    for i=1:length(uci)
        plot(i,meanExpr(cidx==uci(i)),...
            'Marker','o','MarkerSize',4,'markerfacecolor','k','markeredgecolor','none');
    end
    plot(uci,m,...
            'Marker','o','MarkerSize',8,'markerfacecolor','r','markeredgecolor','none','linestyle','none');
    plot(uci,repmat(muData,size(uci)),'--');
    hold off
    ylabel('\mu_{\mu} of cluster');
    xlabel('cluster #');
    xlim([0 length(uci)+1]);
end
% % rescale standard deviation
% term1 = s/max(s)*100;
term1 = (s-sigData)/max(s-sigData);
% term2 = m.*term1;
% nre = clustObj.numObjInClust./sum(clustObj.numObjInClust);
term2 = (m-muData); %term2(term2 < 0) = 0;
% term3 = term1./term2;
cutOff = term1 - term2;
% cutOff = term1;
cutOff = (cutOff-min(cutOff))/max(cutOff-min(cutOff))*100;
cx = find(cutOff==100);
[c,ic] = sort(cutOff,'ascend');
uci = uci(ic);

v = clustObj.Data;
for i=1:length(uci)
    gid = clustObj.cindex==i;
    temp = v(gid,:);
    temp(temp==-inf) = [];
    temp = reshape(temp,numel(temp),1);
    if i==cx
        thrval(i)=muData;
    else
        thrval(i) = prctile(temp,cutOff(i));
    end
end

if figFlag
    subplot(2,2,3)
    plot(1:1:length(uci),c,'ob');
    xlim([0 length(uci)+1]);
    xticks(1:1:length(uci));
    xticklabels(uci);
    ylabel('threshold for cluster');
    xlabel('cluster #');
    
    figure;
    hold on
    yyaxis left
    for i=1:length(uci)
        plot(i,stdExpr(cidx==uci(i))./meanExpr(cidx==uci(i)),...
            'Marker','o','MarkerSize',5,'markerfacecolor','k','markeredgecolor','none');
    end
    ylabel('\sigma/\mu of genes in cluster');
    yyaxis right
    plot(uci,mu_sm,...
            'Marker','o','MarkerSize',10,'markerfacecolor','r','markeredgecolor','none','linestyle','none');
    hold off
    ylabel('\sigma/\mu of cluster');
    xlabel('cluster #');
    xlim([0 length(uci)+1]);
    
    
    figure;
    % figure out number of number of columns and rows in figure
    nrow = 2;
    ncol = 2;
    while nrow*ncol <= length(uci)
        ncol = ncol + 1;
        if nrow*ncol <= length(uci)
            nrow = nrow + 1;
        end
    end
    v = clustObj.Data;
    geneTisMat = zeros(length(uci),(length(xedge)-1));
    for i=1:length(uci)
        gid = clustObj.cindex==i;
        temp = v(gid,:);
        temp(temp==-inf) = [];
        temp = reshape(temp,numel(temp),1);
        geneTisMat(i,:) = histcounts(temp,xedge);
        subplot(nrow,ncol,i);
        hold on
        histogram(temp,xedge,'edgecolor','none');
%         thrval(i) = prctile(temp,cutOff(i));
        yl = ylim;
        plot([thrval(i) thrval(i)],yl); hold off
        ax = gca; ax.TickDir = 'out';
        xlim([min(xedge)-0.5 max(xedge)+0.5]); xticks(xedge);
        ylabel('# of tissue-genes');
        xlabel('Bins of log_{10}[expr.]');
        set(gca,'fontsize',8);
        title(['Cluster # ' num2str(i)]);
    end
end
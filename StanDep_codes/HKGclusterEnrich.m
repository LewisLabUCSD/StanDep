function [pEnr,cnts] = HKGclusterEnrich(clustObj,hkGenes)

% USAGE:
% % [pEnr,cnts] = HKGclusterEnrich(clustObj,hkGenes)

% INPUTS:
% % clustObj:   cluster objected created by geneExprDist_hierarchy
% % hkGenes:    list of housekeeping genes

% OUTPUTS:
% % pEnr:       pvalues for each cluster indicating whether housekeeping
                % % genes were enriched in that cluster
% % cnts:       coverage of housekeeping genes in each cluster

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

nc = size(clustObj.C,1);
ia = ismember(clustObj.objects,hkGenes);
cidxHKG = clustObj.cindex(ia);
[cnts,x] = histcounts(cidxHKG,1:1:nc+1);

% % perform hypergeometric test
pEnr = ones(nc,1);
% number of house-keeping genes (how many can be hits in the box?)
K = sum(ia);
fprintf('No. of housekeeping genes found in the data = %d\n',K);
% total number of genes (how many are in the box?)
M = length(clustObj.objects);
fprintf('List of enriched clusters:\n');
for j=1:nc
    % number of genes in that cluster (how many we pick?)
    N(j) = sum(clustObj.cindex==j);
    % number of house-keeping genes in that cluster (how many are hits?)
    x = cnts(j);
    % calculate hypergeometric p-value
    pEnr(j) = 1 - hygecdf(x,M,K,N(j));
    if pEnr(j) < 0.05
        fprintf('Cluster # %d\n',j);
    end
end
nc_lab = 1:1:nc;
cnts_enr = cnts; cnts_enr(pEnr > 0.05)=0;
fig = figure;
set(fig,'defaultAxesColorOrder',[0 0 1;1 0 0]);
yyaxis left
hold on
bar(1:1:nc,cnts,'barwidth',1,'edgecolor','none','facecolor','k');
bar(1:1:nc,cnts_enr,'barwidth',1,'edgecolor','none','facecolor','b');
hold off
ax = gca; ax.TickDir = 'out';
ylabel('# of house keeping genes');
yyaxis right
plot(1:1:nc,cnts./N,'marker','o','linestyle','none','markerfacecolor','r','markersize',8);
ax = gca; ax.TickDir = 'out';
xlim([0+0.5 nc+0.5]);
xticks(1:1:nc);
xticklabels(nc_lab);
xlabel('cluster #');
ylabel('fraction of cluster-size are house keeping genes');
grid on
title('clusters enriched in house keeping genes (HYPERGEO p-value < 0.05)');
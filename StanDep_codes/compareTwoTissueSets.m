function [J,J2d,J1,J2] = compareTwoTissueSets(set1,set2,typeFlag,tissues,figFlag)

% USAGE:
% % [J,J2d,J1,J2] = compareTwoTissueSets(set1,set2,typeFlag,tissues,figFlag)

% INPUTS:
% % set1:       inputs for set1 
% % set2:       inputs for set2
% % typeFlag:   format of input sets
                % % list, if set1/2 is a cell array of cell arrays.
                % % matrix, if set1/2 is a double(binary) or logical matrix
                % % struct-rxns', if set1/2 is a cell array of COBRA models and reaction similarity is desired
                % % struct-rxns', if set1/2 is a cell array of COBRA models and gene similarity is desired
% % figFlag:    true, plot the figure, else false

% OUTPUTS:
% % J:          Jaccard similarity matrix across sets for same condition
% % J2d:        Jaccard similarity matrix across sets and conditions
% % J1:         Jaccard similairt matrix across conditions for set1
% % J2:         Jaccard similairt matrix across conditions for set2

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

tissues = strrep(tissues,'_','-');
if strcmp(typeFlag,'list')
    J = zeros(length(set1),1);
    J2d = zeros(length(set1),length(set2));
elseif strcmp(typeFlag,'matrix') || strcmp(typeFlag,'struct-rxns') || strcmp(typeFlag,'struct-genes')
    J = zeros(size(set1,2),1);
    J2d = zeros(size(set1,2),size(set2,2));
end
if strcmp(typeFlag,'list')
    for i=1:length(set1)
        J(i) = length(intersect(set1{i},set2{i}))/length(union(set1{i},set2{i}));
        for j=1:length(set2)
            J2d(i,j) = length(intersect(set1{i},set2{j}))/length(union(set1{i},set2{j}));
        end
    end
elseif strcmp(typeFlag,'matrix')
    for i=1:size(set1,2)
        J(i) = sum(set1(:,i) & set2(:,i))/sum(set1(:,i) | set2(:,i));
        for j=1:size(set2,2)
            J2d(i,j) = sum(set1(:,i) & set2(:,j))/sum(set1(:,i) | set2(:,j));
        end
    end
elseif strcmp(typeFlag,'struct-rxns')
    for i=1:size(set1,2)
        J(i) = length(intersect(set1{i}.rxns,set2{i}.rxns))/length(union(set1{i}.rxns,set2{i}.rxns));
        for j=1:size(set2,2)
            J2d(i,j) = length(intersect(set1{i}.rxns,set2{j}.rxns))/length(union(set1{i}.rxns,set2{j}.rxns));
        end
    end
elseif strcmp(typeFlag,'struct-genes')
    for i=1:size(set1,2)
        J(i) = length(intersect(set1{i}.genes,set2{i}.genes))/length(union(set1{i}.genes,set2{i}.genes));
        for j=1:size(set2,2)
            J2d(i,j) = length(intersect(set1{i}.genes,set2{j}.genes))/length(union(set1{i}.genes,set2{j}.genes));
        end
    end
end
J1 = calcJaccardSimilarity(set1,tissues,typeFlag,false);
J2 = calcJaccardSimilarity(set2,tissues,typeFlag,false);
if figFlag
    figure;
    bar(J);
    xlim([0.3 size(set1,2) + 0.7]); xticks(1:1:size(set1,2)); xticklabels(tissues); xtickangle(45);
    xlabel('cell line #'); ylabel(['Jaccard similarity, \mu_J = ' num2str(mean(J))]);
    figure;
    Z1 = linkage(J1,'complete');
    Y1 = pdist(J1,'euclidean');
    leafOrderg = optimalleaforder(Z1,Y1);
    fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
        'complete','euclidean',cophenet(Z1,Y1));
    subplot(1,2,1);
    [~,~,a] = dendrogram(Z1,0,'Orientation','left','labels',tissues,'Reorder',leafOrderg);
    
    Z2 = linkage(J2,'complete');
    Y2 = pdist(J2,'euclidean');
    leafOrderg = optimalleaforder(Z2,Y2);
    fprintf('Cophenetic correlation coeffcient using %s linkage and %s distance = %0.4f\n',...
        'complete','euclidean',cophenet(Z2,Y2));
    subplot(1,2,2);
    [~,~,b] = dendrogram(Z2,0,'Orientation','left','labels',tissues,'Reorder',leafOrderg);
    figure;
    imagesc(J2d(a,b)); hold on
    [~,ia,ib] = intersect(tissues(a),tissues(b));
    [ib,ic] = sort(ib,'ascend'); ia = ia(ic);
    plot(ib,ia,'Marker','o','MarkerSize',8,'Color','k','LineStyle','none'); hold off
    colorbar;
    xlabel('SET 2');
    ylabel('SET 1');
    ax = gca; ax.TickDir = 'out'; ax.DataAspectRatioMode = 'manual';
    xticks(1:1:length(tissues)); yticks(1:1:length(tissues));
    xticklabels(tissues(b));
    yticklabels(tissues(a));
    xtickangle(45);
end
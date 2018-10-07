function [enzymeData] = comparePromiscuousSpecific(spec,prom,modelData,varargin)

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

if to_plot
    figure;
else
    figure('visible','off');
end
hold on
histogram(log10(specExprVec),'BinEdges',[-8:0.25:5],'DisplayStyle','stairs');
histogram(log10(promExprVec),'BinEdges',[-8:0.25:5],'DisplayStyle','stairs');
hold off
legend('Specialist Enzymes','Promiscuous Enzymes','Location','NorthWest');

% % create a matrix of [spec-prom]-by-prom enzyme connectivity
% specPromMat = zeros(length(spec.enzymes)+length(prom.enzymes),length(prom.enzymes));
% specMatrix = zeros(length(spec.enzymes),length(prom.enzymes));
% for i=1:length(prom.enzymes)
%     rxns = prom.rxns{i};
%     if sum(ismember(spec.rxns,rxns))~=0
%         specPromMat(1:length(spec.enzymes),i) = ismember(spec.rxns,rxns);
%         specMatrix(:,i) = ismember(spec.rxns,rxns);
%     end
% end
% specOnlyExpr = specExprMatrix(sum(specMatrix,2)==0,:);
% temp = specOnlyExpr;
% temp(temp==0) = [];
% temp(temp==-1) = [];
% temp = log10(temp);
% figure;
% histogram(temp,'BinEdges',-6:0.25:6,'DisplayStyle','stairs',...
%     'EdgeColor','b','FaceAlpha',1); hold on
% specExpr = specExprMatrix(sum(specMatrix,2)~=0,:);
% temp = specExpr;
% temp(temp==0) = [];
% temp(temp==-1) = [];
% temp = log10(temp);
% histogram(temp,'BinEdges',-6:0.25:6,'DisplayStyle','stairs',...
%     'EdgeColor','r','FaceAlpha',1); hold off
% legend('E_S','E_S.E_P');
% 
% promMatrix = zeros(length(prom.enzymes),length(prom.enzymes));
% for i=1:length(prom.enzymes)-1
%     rxns = prom.rxns{i};
%     for j=i+1:length(prom.enzymes)
%         if sum(ismember(prom.rxns{j},rxns))~=0
%             specPromMat(length(spec.enzymes)+i,j) = sum(ismember(prom.rxns{j},rxns));%2;
%             specPromMat(length(spec.enzymes)+j,i) = sum(ismember(prom.rxns{j},rxns));%2;
%             promMatrix(i,j) = sum(ismember(prom.rxns{j},rxns));
%             promMatrix(j,i) = sum(ismember(prom.rxns{j},rxns));
%         end
%     end
% end
% promOnlyExpr = promExprMatrix(sum(promMatrix,2)==0,:);
% temp = promOnlyExpr;
% temp(temp==0) = [];
% temp(temp==-1) = [];
% temp = log10(temp);
% figure;
% histogram(temp,'BinEdges',-6:0.25:6,'DisplayStyle','stairs',...
%     'EdgeColor','b','FaceAlpha',1); hold on
% promExpr = promExprMatrix(sum(promMatrix,2)~=0,:);
% temp = promExpr;
% temp(temp==0) = [];
% temp(temp==-1) = [];
% temp = log10(temp);
% histogram(temp,'BinEdges',-6:0.25:6,'DisplayStyle','stairs',...
%     'EdgeColor','r','FaceAlpha',1); hold off
% legend('E_P','E_S.E_P');

enzymeData.enzyme = [spec.enzymes;prom.enzymes];
enzymeData.value = enzAbund;
% enzymeData.specPromMat = specPromMat;
enzymeData.rxns = [spec.rxns;prom.rxns];
enzymeData.Tissue = modelData.Tissue;
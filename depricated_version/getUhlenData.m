function [expressionData] = getUhlenData(geneList)

fileName = 'rna_tissue_fpkm.txt';
T = readtable(fileName,'Delimiter','\t','ReadVariableNames',false);
eval(strcat(regexprep(T.Var1{1,1},'\.| ',''),' = T.Var1(2:end);'));
% for i=1:32
%     eval(strcat(regexprep(T.Var3{1,1},'\.| ',''),' = str2double(T.Var3(2:end));'));
% end
eval(strcat(regexprep(T.Var19{1,1},'\.| ',''),' = str2double(T.Var19(2:end));'));
eval(strcat(regexprep(T.Var20{1,1},'\.| ',''),' = str2double(T.Var20(2:end));'));
eval(strcat(regexprep(T.Var21{1,1},'\.| ',''),' = str2double(T.Var21(2:end));'));
eval(strcat(regexprep(T.Var22{1,1},'\.| ',''),' = str2double(T.Var22(2:end));'));
eval(strcat(regexprep(T.Var23{1,1},'\.| ',''),' = str2double(T.Var23(2:end));'));
eval(strcat(regexprep(T.Var24{1,1},'\.| ',''),' = str2double(T.Var24(2:end));'));
eval(strcat(regexprep(T.Var25{1,1},'\.| ',''),' = str2double(T.Var25(2:end));'));
eval(strcat(regexprep(T.Var26{1,1},'\.| ',''),' = str2double(T.Var26(2:end));'));
eval(strcat(regexprep(T.Var27{1,1},'\.| ',''),' = str2double(T.Var27(2:end));'));
eval(strcat(regexprep(T.Var28{1,1},'\.| ',''),' = str2double(T.Var28(2:end));'));
eval(strcat(regexprep(T.Var29{1,1},'\.| ',''),' = str2double(T.Var29(2:end));'));
eval(strcat(regexprep(T.Var30{1,1},'\.| ',''),' = str2double(T.Var30(2:end));'));
eval(strcat(regexprep(T.Var31{1,1},'\.| ',''),' = str2double(T.Var31(2:end));'));
eval(strcat(regexprep(T.Var32{1,1},'\.| ',''),' = str2double(T.Var32(2:end));'));
eval(strcat(regexprep(T.Var33{1,1},'\.| ',''),' = str2double(T.Var33(2:end));'));
eval(strcat(regexprep(T.Var34{1,1},'\.| ',''),' = str2double(T.Var34(2:end));'));
eval(strcat(regexprep(T.Var35{1,1},'\.| ',''),' = str2double(T.Var35(2:end));'));
eval(strcat(regexprep(T.Var36{1,1},'\.| ',''),' = str2double(T.Var36(2:end));'));
eval(strcat(regexprep(T.Var37{1,1},'\.| ',''),' = str2double(T.Var37(2:end));'));
eval(strcat(regexprep(T.Var38{1,1},'\.| ',''),' = str2double(T.Var38(2:end));'));
eval(strcat(regexprep(T.Var39{1,1},'\.| ',''),' = str2double(T.Var39(2:end));'));
eval(strcat(regexprep(T.Var40{1,1},'\.| ',''),' = str2double(T.Var40(2:end));'));
eval(strcat(regexprep(T.Var41{1,1},'\.| ',''),' = str2double(T.Var41(2:end));'));
eval(strcat(regexprep(T.Var42{1,1},'\.| ',''),' = str2double(T.Var42(2:end));'));
eval(strcat(regexprep(T.Var43{1,1},'\.| ',''),' = str2double(T.Var43(2:end));'));
eval(strcat(regexprep(T.Var44{1,1},'\.| ',''),' = str2double(T.Var44(2:end));'));
eval(strcat(regexprep(T.Var45{1,1},'\.| ',''),' = str2double(T.Var45(2:end));'));
eval(strcat(regexprep(T.Var46{1,1},'\.| ',''),' = str2double(T.Var46(2:end));'));
eval(strcat(regexprep(T.Var47{1,1},'\.| ',''),' = str2double(T.Var47(2:end));'));
eval(strcat(regexprep(T.Var48{1,1},'\.| ',''),' = str2double(T.Var48(2:end));'));
eval(strcat(regexprep(T.Var49{1,1},'\.| ',''),' = str2double(T.Var49(2:end));'));
eval(strcat(regexprep(T.Var50{1,1},'\.| ',''),' = str2double(T.Var50(2:end));'));

tissueNames = setdiff(who,{'T';'fileName';'ensg_id';'geneList'});
rawData = eval(strcat('[',strjoin(tissueNames),']'));

if nargin < 1
    geneList = ensg_id;
end
data = zeros(numel(geneList),numel(tissueNames));
if nargin < 1
    data = rawData;
    expressionData.genesymbol = T.Var3(2:end);
else
    c = 0;
    for i=1:length(geneList)
        if sum(ismember(ensg_id,geneList{i}))~=0
            c = c+1;
            g{c,1} = geneList{i};
            data(c,:) = rawData(ismember(ensg_id,geneList{i}),:);
        end
    end
    ensg_id = g;
end
expressionData.gene = ensg_id;
expressionData.Tissue = strrep(tissueNames,'fpkm','');
expressionData.valuebyTissue = data;
expressionData.unit = 'FPKM';
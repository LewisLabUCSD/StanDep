function [rnaData,tissues] = getNCI60Data

filename = 'nci60_RNA__Agilent_mRNA_log2.txt';
opts = detectImportOptions(filename);
T = readtable(filename,opts);
rnaData.Tissue = opts.VariableNames(3:end)';
tissues = rnaData.Tissue;
rnaData.gene = T.GeneId;
for i=1:length(rnaData.Tissue)
    rnaData.valuebyTissue(:,i) = eval(strcat('2.^T.',rnaData.Tissue{i},';'));
end
startPos = regexp(rnaData.Tissue,'_','once');
rnaData.Tissue = cellfun(@extractAfter,rnaData.Tissue,startPos,'uniformoutput',false);
rnaData.Tissue = linearization_index(rnaData.Tissue,'cols');
rnaData.Tissue = strrep(rnaData.Tissue,'_','');
gene = rnaData.gene;
gene = num2str(gene);
gene = cellstr(gene);
gene = strrep(gene,' ','');
rnaData.gene = gene;

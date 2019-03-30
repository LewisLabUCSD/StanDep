function [nRxns,entropy_score] = simulation_analysis(sim_output,selection_mat,rnaData,model_file,enzymeData,f_save)

load(model_file,'model');
%% Calculate which reactions are not associated with a gene 
rxns = enzymeData.rxns(sum(enzymeData.value,2)~=0);
rxns = unique([linearization_index(rxns(cellfun(@iscell,rxns)),'cols'); rxns(~cellfun(@iscell,rxns))]);
ix = ismember(model.rxns,rxns); %ix is are rxns associated with a gene

%% Calculate a. Jaccard Index for each sample b. Number of rxns kept for each sample
num_samples = size(sim_output,2);

entropy_score = cell([num_samples,1]); %zeros(num_samples); 
nRxns = cell([1,num_samples]); %zeros(num_samples,nSim);
for i = 1:num_samples
    % For each tissue from sim_output only the values where nSim is 1
    sim_inds = selection_mat(:,i) == 1;
    curr_sample = sim_output(:,i,sim_inds);
    curr_sample = reshape(curr_sample, [size(sim_output,1),sum(sim_inds)]);
    
    % Remove reactions nont associated with a gene
    curr_sample = curr_sample(ix,:);
    
    % Total number of reactions kept 
    nRxns{i} = sum(curr_sample);
    
    % Jaccard correlation
    entropy_score{i} = calc_entropy(curr_sample);
end

%% Get the actual values of number of reactions
[data_mat] = funcConstructModels(rnaData,model_file); % Get original value
data_mat = data_mat(ix,:);
actual_vals = zeros([size(data_mat,2),1]);
for i = 1:size(data_mat,2)
   actual_vals(i) = sum(data_mat(:,i));
end



%% Number of reactions for each tissue boxplot
box_vector = cell2mat(nRxns);
box_vector = transpose(box_vector);
name_vector = strings(size(box_vector));
curr_count = 1;
xticklabel_names = cell([1,num_samples]);
for i =1:num_samples
    name_vector(curr_count: curr_count + length(nRxns{i})-1) = i; % rnaData.Tissue{i};
    xticklabel_names{i} = rnaData.Tissue{i}; % For labeling
    curr_count = curr_count + length(nRxns{i});
end
name_vector = char(name_vector); % Padding the strings to have same length

% Plot
f = figure;
boxplot(box_vector,name_vector);
title('Number of reactions kept across all simulations')
xlabel('Tissue')
ylabel('Number of Reactions')

xt = xticks;
xticklabels(xticklabel_names)
xtickangle(90)

%Plot the actual values
hold on 
scatter(xt, actual_vals,20,'filled','green')
savefig(strcat(f_save,'_numRxns'))
saveas(gcf,strcat(f_save,'_numRxns'),'png')

%% Entropy score boxplot
entropy_vector = cell2mat(entropy_score);
%entropy_vector = transpose(entropy_vector);
name_vector = strings(size(entropy_vector));
curr_count = 1;
for i =1:num_samples
    name_vector(curr_count: curr_count + length(nRxns{i})-1) = rnaData.Tissue{i};
    curr_count = curr_count + length(nRxns{i});
end
name_vector = char(name_vector); % Padding the strings to have same length

% Plot
figure;
boxplot(entropy_vector,name_vector);
title('Entropy of each reaction across all simulations')
xlabel('Tissue')
ylabel('Entropy Score')
xtickangle(90)
savefig(strcat(f_save,'_entropyBoxplot'))
saveas(gcf,strcat(f_save,'_entropyBoxplot'),'png')

%% Entropy histogram for each tissue 
%p = 0.9  entropy = 0.469

% determine shape of the figure- create a square
nrows = 2;
ncols = 2;
while nrows*ncols <= num_samples
    ncols = ncols + 1;
    if nrows*ncols <= num_samples
        nrows = nrows + 1;
    end
end

figure
for i = 1:num_samples
    subplot(nrows,ncols,i)
    histogram(entropy_score{i})
    %[N,edges] = histcounts(entropy_score{i});
    %bar(log10(N));
    set(gca,'YScale','log')
    xlabel('Entropy score')
    ylabel('Frequency')
    title(rnaData.Tissue{i})
end
savefig(strcat(f_save,'_entropyHist'))
saveas(gcf,strcat(f_save,'_entropyHist'),'png')

end
function [entropy_score] = calc_entropy(dat)
%% Computes entropy for binary matrix: 
%  calculates row-wise the entropy and then takes the average
p = sum(dat,2)/size(dat,2); %sum(dat)/length(dat);
q = 1-p;
entropy_score = -(log2(p).*p + log2(q).*q);
entropy_score(isnan(entropy_score)) = 0;
%entropy_score = mean(entropy_score);
end
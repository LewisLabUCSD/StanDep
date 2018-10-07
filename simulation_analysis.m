function [nRxns,entropy_score] = simulation_analysis(sim_output,selection_mat,rnaData)

%% Calculate a. Jaccard Index for each sample b. Number of rxns kept for each sample
num_samples = size(sim_output,2);

entropy_score = cell([num_samples,1]); %zeros(num_samples); 
nRxns = cell([1,num_samples]); %zeros(num_samples,nSim);
for i = 1:num_samples
    % For each tissue from sim_output only the values where nSim is 1
    sim_inds = selection_mat(:,i) == 1;
    curr_sample = sim_output(:,i,sim_inds);
    curr_sample = reshape(curr_sample, [size(sim_output,1),sum(sim_inds)]);
    
    % Total number of reactions kept 
    nRxns{i} = sum(curr_sample);
    
    % Jaccard correlation
    entropy_score{i} = calc_entropy(curr_sample);
end

%% Boxplots for each tissue
box_vector = cell2mat(nRxns);
box_vector = transpose(box_vector);
name_vector = strings(size(box_vector));
curr_count = 1;
for i =1:num_samples
    name_vector(curr_count: curr_count + length(nRxns{i})-1) = rnaData.Tissue{i};
    curr_count = curr_count + length(nRxns{i});
end

name_vector = char(name_vector); % Padding the strings to have same length
figure;
boxplot(box_vector,name_vector);
title('Number of reactions kept across all simulations')
xlabel('Tissue')
ylabel('Number of Reactions')

%% Boxplots for each entropy score
entropy_vector = cell2mat(entropy_score);
%entropy_vector = transpose(entropy_vector);
name_vector = strings(size(entropy_vector));
curr_count = 1;
for i =1:num_samples
    name_vector(curr_count: curr_count + length(nRxns{i})-1) = rnaData.Tissue{i};
    curr_count = curr_count + length(nRxns{i});
end

name_vector = char(name_vector); % Padding the strings to have same length
figure;
boxplot(entropy_vector,name_vector);
title('Entropy of each reaction across all simulations')
xlabel('Tissue')
ylabel('Entropy Score')
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
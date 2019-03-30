function [sim_output,selection_mat] = subset_simulation(rnaData,model_file, varargin)
% Function that will subsample the columns of a d-by-nSamples 2D array,run a
% function and compute the outputs correlation across all subsamples
% 
% Optional variables: 
% nSim: Number of simulations to runn
% fract_samples: Fraction of the samples to subsample
% save_f: Output file name to save variables. If not put, will not save.

numvarargs = length(varargin);
if numvarargs > 3
    error('TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {10,0.5,''};
optargs(1:numvarargs) = varargin;
% Place optional args as variable names
[nSim,frac_samples,save_f] = optargs{:};

%% Run constructModel over entire samples first
% This will tell you how many dimensions there are
[data_mat] = funcConstructModels(rnaData,model_file);
close all
%% Initialize variables
% Initialize selection_mat, an nSim-by-nSamples 2D array which is an indicator matrix for
% which samples will be selected in each simulation
num_samples = size(data_mat,2);
selection_mat = zeros(nSim,num_samples);
subsample_num = floor(frac_samples*num_samples);
disp(['Number of samples for each simulation: ',num2str(subsample_num)]);

% Initialize sim_output, a nRxns-by-nSamples-by-nSim 3D array
sim_output = -1*ones([size(data_mat),nSim]); %Binary output but make samples that were never seenn -1


parfor i =1:nSim 
    disp(['Simulation Number: ', num2str(i)]);
    select_inds = randsample(num_samples,subsample_num);
    
    %%% Variable assignment selectionn_mat
    v = zeros(1, num_samples); %This is for paallelization (PCT)
    v(:,select_inds) = 1;
    selection_mat(i,:) = v; 
    %selection_mat(i,select_inds) = 1; %CANT USE this in PCT
    %%% 
    
    currRNA = rnaData; 
    currRNA.Tissue = currRNA.Tissue(select_inds);
    currRNA.valuebyTissue = currRNA.valuebyTissue(:,select_inds);
    % Run function on currRNA and store output in sim_output
    curr_out = funcConstructModels(currRNA, model_file); %func(curr_mat);
%     disp('curr out size')
%     disp(size(curr_out));
    close all
    
    %%% Variable assignment sim_output
    w = -1*ones(size(data_mat));    % For PCT
    w(:,select_inds) = curr_out;
    sim_output(:,:,i) = w
    %sim_output(:,select_inds,i) = curr_out; %CANT USE this in PCT
    %%% 
    
end

if ~strcmp(save_f,'')
    save(save_f,'sim_output','selection_mat','nSim','frac_samples')
end

% %% Calculate Jaccard Index for each sample
% jaccard = zeros(num_samples); 
% for samp = 1:num_samples
%     % For each tissue from sim_output only the values where nSim is 1
%     samp_inds = selection_mat(:,samp) == 1;
%     curr_sample = sim_output(:,samp,samp_inds);
% end

end

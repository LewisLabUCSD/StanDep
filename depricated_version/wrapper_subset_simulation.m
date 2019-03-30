%Wrapper
rnaData = getUhlenData;
model_file = 'ENSG_Recon22_Mapping';
%load('ENSG_Recon22_Mapping','model');

if exist('p','var')
    delete(p);
end
    
p = parpool(4);

tic
nSim = 4;% 1000;
frac_samples = 0.75;
save_f = 'subset_results/subset_simulation_Uhlen_003.mat';
[sim_output,selection_mat] = subset_simulation(rnaData,model_file, nSim, frac_samples,save_f);

toc; 

delete(p)


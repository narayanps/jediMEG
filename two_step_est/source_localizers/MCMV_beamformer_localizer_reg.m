function [loc_est, id_est, G_mcmv, data_cov_inv] = MCMV_beamformer_localizer_reg(Y, noise_cov, G, pnts, included_sp, num_sources)
%addpath('/m/nbe/scratch/braintrack/pnas_mne_results')
num_trials = size(Y,3);
T = size(Y,2);
data_cov = compute_cov(Y, num_trials, T);


%%%%%%%%%%%reg%%%%%%%%%%%%%%%%%%%%%
reg=0.05;
alpha=reg*trace(data_cov)/length(data_cov);
nchan = size(data_cov,2);
data_cov_reg=data_cov+alpha*eye(nchan);

alpha=reg*trace(noise_cov)/length(noise_cov);
nchan = size(noise_cov,2);
noise_cov_reg=noise_cov+alpha*eye(nchan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:size(G,2)
    eta_ai(i) =  compute_ai(G(:,i), noise_cov_reg , data_cov_reg);
end

loc(1) = find(eta_ai == max(eta_ai));

for j=2:1:num_sources

    for i=1:1:size(G,2)
        id = find(loc == i);
        if isempty(id)
            eta_mai(i) = compute_mai([G(:,loc) G(:,i)], noise_cov_reg , data_cov_reg);
        end
    end

eta_mai(loc) = 0;
loc(j) = find(eta_mai == max(eta_mai));
clear eta_mai
end

loc_est = reshape(pnts(included_sp(loc),:)', 3*num_sources,1);
G_mcmv = G(:,loc);
id_est = loc;
data_cov_inv = inv(data_cov_reg);
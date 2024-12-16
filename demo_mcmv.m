function demo_mcmv(sim,type)
path_ = pwd;
addpath(genpath(strcat(path_,'/MVAR_simu')));
addpath(genpath(strcat(path_,'/two_step_est')));
addpath(genpath(strcat(path_,'/misc')))

meas_snr=5;
bio_snr_vals=[1 3 5 10];
bio_snr_arr=repmat(bio_snr_vals, [50 1]);
T=5000;
num_bio_ns=500;

num_sources=3;
dist_thr = 30;

%load head model
load simu_data.mat
[source_inds] = select_source_inds(pnts, num_sources, dist_thr);
[q, A] = three_source_model(T, type);
[Y, sigma_b, sigma_m, est_snr] = gen_MEG_simu(q,G,source_inds, meas_snr, bio_snr_arr(sim), num_bio_ns);

% ground truth
GT.A=A;
GT.source_inds=source_inds;
GT.sigma_m=sigma_m;
GT.q=q;
GT.Y=Y;
GT.source_loc = pnts(source_inds,:);

%prepare model struct for estimation
load estim_data.mat
model.G=G_e;
num_sources=3;

[loc_est,id_est,G_mcmv, data_cov_inv] = MCMV_beamformer_localizer_reg(Y,eye(size(Y,1)), model.G, pnts_e, 1:length(pnts_e), num_sources);

for tr=1:1:size(Y,3)
 wT = inv(G_mcmv(:, :)' * data_cov_inv * G_mcmv(:, :)) * G_mcmv(:, :)' * data_cov_inv;
 amp_est(:,:,tr) = wT*squeeze(Y(:,:,tr));
end
MCMV_est.loc_est=loc_est;
MCMV_est.amp_est = amp_est;

end





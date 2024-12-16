function demo_lcmv(sim,type)
path_ = pwd;
addpath(genpath(strcat(path_,'/MVAR_simu')));
addpath(genpath(strcat(path_,'/two_step_est')));
addpath(genpath(strcat(path_,'/misc')))

meas_snr=5;
bio_snr_vals=[1 3 5 10];
bio_snr_arr=repmat(bio_snr_vals, [50 1]);
T=5000;
num_bio_ns=1000;
num_sources = 3;
dist_thr = 30;


%load head model for simu - lead-filed G and mesh points pnts for simulation
load simu_data

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

%prepare model struct for estimation - lead-field and mesh points onts_e for estimation
load estim_data

num_sources=3;
model.G = G_e; %different lead-field for estimation

num_trials = size(Y,3);
T = size(Y,2);
C = compute_cov(Y, num_trials, T);

for i=1:size(model.G,2)
    [eta_ai(i),Creg] =  compute_LCMV_reg(model.G(:,i), 1.*eye(size(Y,1)) , C);
end

[~, id2] = sort(eta_ai, 'descend');
ids=id2(1:num_sources);
loc_est=pnts_e(ids,:);
G_nai=model.G(:,ids);
Cinv = inv(Creg);

for tr=1:1:size(Y,3)
  wT = inv(G_nai(:, :)' * Cinv * G_nai(:, :)) * G_nai(:, :)' * Cinv;
  amp_est(:,:,tr) = wT*squeeze(Y(:,:,tr));
end

LCMV_est.loc_est=loc_est;
LCMV_est.amp_est = amp_est;


end





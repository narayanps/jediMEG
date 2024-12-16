function [ind_, r_, G_] = set_part_traj_mcmv_prior(Y, noise_cov, G, pnts_e, included_sp, num_sources, T)
[loc_est, id_est, G_mcmv, ~] = MCMV_beamformer_localizer(Y, noise_cov, G, pnts_e, included_sp, num_sources);

for j=1:1:num_sources
    ind_(j,:) = repmat(id_est(j), [1 T]);
end
r_ = repmat(loc_est,[1 T]);
G_ =repmat(G_mcmv,[1 1 T]);


end

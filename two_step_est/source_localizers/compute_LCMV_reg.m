function [eta, data_cov_reg] = compute_LCMV_reg(h, noise_cov, data_cov)
reg=0.05;
alpha=reg*trace(data_cov)/length(data_cov);
nchan = size(data_cov,2);
data_cov_reg=data_cov+alpha*eye(nchan);

alpha=reg*trace(noise_cov)/length(noise_cov);
nchan = size(noise_cov,2);
noise_cov_reg=noise_cov+alpha*eye(nchan);

eta = (h' * inv(noise_cov_reg)* h )/ (h' * inv(data_cov_reg)* h );
end

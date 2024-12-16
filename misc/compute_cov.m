function data_cov = compute_cov(Y, num_trials, T)
data_cov=0;
for j=1:1:num_trials
    data_cov = data_cov + (squeeze(Y(:,:,j))*squeeze(Y(:,:,j))');
end
data_cov = data_cov./T;
data_cov=double(data_cov);
end


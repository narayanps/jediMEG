function [id, r_, G_] = set_all_part_traj_v2(pnts, G, dist_thr, num_sources, N, m0, P0)
M=length(pnts);
dist_sources=zeros(1,N);
r_=zeros(3*num_sources,N);
G_=zeros(size(G,1),num_sources,N);
min_dist = 0;
r_1 = repmat(m0, [1 N]) + chol(P0)*randn(3*num_sources, N);
for j=1:N
    pos = reshape(r_1(:,j), [3 num_sources])';
    for i=1:1:size(pos,1)
        [dist_vec, id(i,j)] = find_nearest_id(pnts,pos(i,:));
    end
G_(:,:,j) = G(:, id(:,j)');
r_(:,j) = reshape(pnts(id(:,j),:)', [3*num_sources,1]);
end

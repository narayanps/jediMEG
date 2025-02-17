function [q, P, Kf] = gpu_kalman_update(q, P, Y_, G, G_, XI, N, Jt)
gputimes = @(A, B) pagefun(@mtimes, A, B);      % A*B on GPU
gpurdivide = @(A, B) pagefun(@mrdivide, A, B);  % A/B on GPU
gputranspose = @(A) pagefun(@transpose, A);     % A'
ny=size(Y_,1);
nq=size(q,1);
ns=size(G,2);
Kf = gpurdivide(gputimes(P(:,1:ns,:), gputranspose(G)), gputimes(gputimes(G,P(1:ns,1:ns,:)), gputranspose(G)) + XI);
%Kf = gpurdivide(gputimes(P, gputranspose(G)), gputimes(gputimes(G,P), gputranspose(G)) + XI);
Kf_ = repmat(Kf, [1 1 1 Jt]);
Kf_=reshape(Kf_, nq,ny,N*Jt);
q_(:,1,:) = q;
v(:,1,:) = Y_ - squeeze(gputimes(G_, q_(1:ns,:,:)));
q = q + squeeze(gputimes(Kf_, v)); 
P = P - gputimes(gputimes(Kf, G), P(1:ns,:,:));
%P = P - gputimes(gputimes(Kf, G), P);
q=gather(squeeze(q));
P=gather(P);
Kf=gather(Kf);
end
            
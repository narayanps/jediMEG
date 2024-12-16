function [state, params, LL_complete] = demo_saem_meg(Y, model)

% Y is the MEG data (real or simulated from ECoG)

% model struct contains
    % G - lead field matrix
    % included_verts - included source points
    % pnts_e - mesh points xyz in mm

path_ = pwd;
addpath(genpath(strcat(path_,'/joint_est')));


%scaling the lead-field to match source activity scale
model.G=double(G)*1e-9;

%SAEM iteration parameters
num_iter = 1000;
b=199;
kappa = 1; 
gamma = zeros(1,num_iter);
gamma(1:2) = 1;
gamma(3:b) = 0.98;
gamma(b+1:end) = 0.98*(((0:num_iter-(b+1))+kappa)/kappa).^(-0.7);

%set number of sources, particles, order etc
Ns=5;
Np=300;
P=14; %model order
whiten_flag=0;
dist_thr = 10; %in mm, min distance between source dipoles in each particle


%initial value for A and V
mean_A0 =zeros(1,Ns)-0.9;
Sigma_A0 = 1e-3*eye(Ns);
lambdamax=10;
while lambdamax >1       
for jj = 1:P
    A_p = mean_A0' + chol(Sigma_A0)*randn(Ns,1);         %diagonal entries of A matrix are U[a,b]
    A_p = diag(diag(A_p));
    A0(Ns*(jj-1)+1:Ns*(jj-1)+Ns,1:Ns) = A_p.*eye(Ns);
end 
A0 = A0';
A0(Ns+1:Ns*P, 1:Ns*(P-1)) = eye(Ns*(P-1));
A0(Ns+1:Ns*P, 1+Ns*(P-1):end) = zeros(Ns*(P-1), Ns);
lambda=eig(squeeze(A0));lambdamax=max(abs(lambda));
end  

%V_init_range = [0.1 0.9];
%V_q0 = diag(V_init_range(1)+(V_init_range(2)-V_init_range(1))*rand(Ns*P,1));
V_q0=eye(Ns*P)  ;
V_q0(Ns+1:Ns*P, 1:Ns)=0;
V_q0(1:Ns*P, Ns+1:end)=0;

% init struct for SAEM algorithm
init.A0 =A0;
init.V0 = V_q0;
init.q0 = zeros(Ns*P,1);
init.P0=1*eye(Ns*P);
init.sigma_m0 = 1; %eye(size(G,1));
init.sigma_b0 = 1;


%opt_params struct
opt_params.Niter = num_iter;
opt_params.gamma = gamma;
Y_avg = squeeze(mean(Y,3)); 
gpu_flag=1;
mcmv_prior_flag=1;

[state, params, LL_complete]= saem(init, model, opt_params, Y_avg, Y, Y, Ns, Np, P, dist_thr, gpu_flag, whiten_flag, mcmv_prior_flag);

end





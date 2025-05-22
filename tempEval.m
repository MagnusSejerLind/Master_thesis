% min(N): RMSE(N)<eps

clc,clear,close all
set(0,'defaultTextInterpreter','latex');



% Number of time steps (dt def. i func.)
l1 = 2:1:99;
l2 = 100:10:300-10;
l3 = 300:100:1000;
N = [l1,l2,l3];


% calc RMSE(N)
for i = 1:length(N)
    disp(N(i))
    [RMSE(i)] = main_rmse_temp(N(i));
end

%%

figure()
plot(N,RMSE,'k.-',LineWidth=2,MarkerSize=20)
xlabel('N')
ylabel('RMSE')
grid
title('RMSE over number of timesteps')



%% function
function [RMSE_tot] = main_rmse_temp(N)
rng('default')

%% System & options
opt.sysType = "chain";  % ["chain" / "frame"] - Type of system
opt.method = "TA";      % ["TA"/"ME"] - Virtuel sensing method (Toeplitz's/Modal expansion)
opt.out_type = 0;       % [0/1/2] - Define output type (0=disp, 1=vel, 2=accel)
opt.numDOF = 6;         % [-int.-] - Number of DOF --ONLY FOR CHAIN SYSTEM
opt.nonlinear = 1;      % [0/1] - Include nonlinearties in the system
opt.nonlinType = 1;     % [0/1] - Define type of nonlineaties (0=constant / 1=varied)
opt.error_mod = 0;      % [0/1] - Include error modeling and noise
opt.psLoads = 1;        % [0/1] - Apply pseodu loads to convert nonlinear system to linear model
opt.condimp = 1;        % [0/1] - Improve Toeplitz's matrix condition by truncation
opt.plot = 0;           % [0/1] - plot results
opt.animate = 0;        % [0/1] - Animates the displacements of the structure
opt.aniSave = 0;        % [0/1] - Save animation
% opt;

in_dof = [1];         % Input DOF
out_dof = [6,2];        % Output DOF
% out_dof = [1 2 3 4 5 6 7 8 9 10 12 15 16 18 19 20 22 23 24];        % Output DOF
% frame dof: (x,y,θ)
%% System modeling

[dof,m,k,xi_int] = systemSetup(opt);
r = numel(in_dof);
ms = numel(out_dof);

% IC
d0 = zeros(dof,1);
v0 = zeros(dof,1);
z0 = [d0;v0];

% Time
% N = 500;  :-function input-:
dt = 0.01;
t = 0:dt:(N-1)*dt;

% Input (dofs defined earlier)
u_mag = 10;
u = ones(r,N)*u_mag;
u = u.*sin(t*5);
% u = zeros(r,N);
% u(N*0.2:N*0.3) = u_mag;
U = u(:);
if opt.error_mod == 1
    % u_acc = u; U_acc = U;
    [~,~,snr] = modeling_error(0,0);
    U = awgn(U,snr,'measured');
    u = reshape(U,r,N);
end


% Actucal system - no mod. error
if opt.sysType == "chain"; [M_acc,~,K_acc] = chain(m,m*0,k,dof); end
addBeamError = 0;
if opt.sysType == "frame"; [M_acc,K_acc,dof,snr] = beamStruc(opt,addBeamError); end
[Phi_acc,Lambda_acc] = eig(K_acc,M_acc);    % modal and spectral matrix
[omegaN_acc,i2] = sort(sqrt(diag(Lambda_acc))); % Natural freq.
omegaN_acc = real(omegaN_acc);
Phi_acc = Phi_acc(:,i2);
dd_acc = sqrt(diag(Phi_acc'*M_acc*Phi_acc));
% aa_acc = Phi_acc*diag(1./dd_acc);    % Mass-normalized Phi (eigenvec.)
[alpha_acc,beta_acc] = raylieghDamp(omegaN_acc,xi_int);
C_acc = alpha_acc*M_acc + beta_acc*K_acc;
% C_modal_acc = round(Phi_acc'*C_acc*Phi_acc,10);


% Base system
if opt.sysType == "chain"
    if opt.error_mod == 1; [k,m,snr] = modeling_error(k,m); end
    [M,~,K] = chain(m,m*0,k,dof);
end
addBeamError = 1;
if opt.sysType == "frame"; [M,K,dof,snr] = beamStruc(opt,addBeamError); end
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi = Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi)); % Mass norm M
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
[alpha,beta] = raylieghDamp(omegaN,xi_int);
C = alpha*M + beta*K;
C_modal = round(Phi'*C*Phi,10);
% xi = diag(C_modal) ./ (2*omegaN);


% Extended system - full output
in_dof_ex = in_dof;
out_dof_ex = (1:1:dof);
dof_ex = numel(out_dof_ex);
r_ex = numel(in_dof_ex);
ms_ex = numel(out_dof_ex);


% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,opt.out_type,dt);
[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof,out_dof_ex,opt.out_type,dt);

% Toeplitz's matricies
[H] = ToeplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);
[H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

% Nonlinearities
if opt.nonlinear == 1
    nl_mag = 0.5;
    if opt.nonlinType == 0
    % varied / const. mag.
        cf_nl = nl_mag;    % coeffcient of nonlinear damping
        kf_nl = nl_mag;    % coeffcient of nonlinear stiffness
    else
        cf_nl = rand(1,dof)*nl_mag;
        kf_nl = rand(1,dof)*nl_mag;
    end
else
    cf_nl = 0;
    kf_nl = 0;
end

%% Compute outputs

% Base system
z_old = z0;
z_new = zeros(size(z_old));
fd_nl = zeros(size(z_old));
fk_nl = zeros(size(z_old));
y = zeros(ms,length(t));
for i = 1:N
    fd_nl(dof+1:end) = cf_nl*z_old(dof+1:end).*abs(z_old(dof+1:end));   % non-linear damping force (velocity dependent)
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);                        % non-linear stiffness force (displacement dependent)

    z_new = Ad*z_old + Bd*u(:,i) - fd_nl - fk_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);
if opt.error_mod == 1
% output noise
    Y = awgn(Y,snr,'measured');
    y = reshape(Y,ms,N);
end

% Actual system
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));
fd_nl_acc = zeros(size(z_old_acc));
fk_nl_acc = zeros(size(z_old_acc));
y_acc = zeros(dof,length(t));
for i = 1:N
    fd_nl_acc(dof+1:end) = cf_nl*z_old_acc(dof+1:end).*abs(z_old_acc(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl_acc(dof+1:end)  = kf_nl*(z_old_acc(1:dof).^3);                            % non-linear stiffness force (displacement dependt)

    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i) - fd_nl_acc - fk_nl_acc;
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);

%% Condition improvement

if opt.condimp == 1
% Removes first ms rows, and last r columns in the Toeplitz's matricies

H = H(ms+1:end, 1:end-r);
H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);

U = U(1:end-r);
Y = Y(ms+1:end);

end

%% Pseudo loads

if opt.psLoads == 1

    % Full input, org. out
    in_dof_FI = 1:1:dof;
    r_FI = numel(in_dof_FI);
    [Ad_FI,Bd_FI,Cd_FI,Dd_FI] = systemMatriciesSS_dis(M,K,C,dof,in_dof_FI,out_dof,opt.out_type,dt);
    [H_FI] = ToeplitzMatrix(N,ms,r_FI,Ad_FI,Bd_FI,Cd_FI,Dd_FI);

    if opt.condimp
        H_FI = H_FI(ms+1:end, 1:end-r_FI);
    end

    % % Tikhonov regularization parameter
    % lambda = 0.071969;                   % adjust this to trade bias vs. variance
    % % lambda = 1e-10;
    % % build regularized normal equations
    % A = H_FI' * H_FI + lambda * eye(size(H_FI,2));
    % b = H_FI' * (Y - H * U);
    % % solve for Gamma
    % Gamma = A \ b;

    Gamma = pinv(H_FI)*(Y - H*U);

else
    Gamma = zeros(dof*length(t),1);
    if opt.condimp == 1 && opt.psLoads ~= 1; Gamma = Gamma(1:end-dof); end
end
%% Output estimation

if opt.method == "TA"
% Toeplitz approach

    % Full input, full output
    in_dof_FIFO = 1:1:dof;
    r_FIFO = numel(in_dof_FIFO);
    ms_FIFO = numel(out_dof_ex);
    [Ad_FIFO,Bd_FIFO,Cd_FIFO,Dd_FIFO] = systemMatriciesSS_dis(M,K,C,dof,in_dof_FIFO,out_dof_ex,opt.out_type,dt);
    [H_FIFO] = ToeplitzMatrix(N,ms_FIFO,r_FIFO,Ad_FIFO,Bd_FIFO,Cd_FIFO,Dd_FIFO);
    % Full output, org. input
    H_ex;

if opt.condimp
    H_FIFO = H_FIFO(ms_FIFO+1:end, 1:end-r_FIFO);
end

    Y_est = H_ex*U + H_FIFO*Gamma;

    if opt.condimp == 1
    y_est = reshape(Y_est, dof, N-1);  % decollapse dof columns
    else
        y_est = reshape(Y_est, dof, N);  % decollapse dof columns
    end
end

if opt.method == "ME"
    % Modal expansion

    mu1 = out_dof;              % Observed nodes 
    mu2 = 1:dof; mu2(mu1)=[];   % Unobserved nodes
    eta1 = 1:numel(mu1);        % Retained modes

    Phi_mu1_eta1 = aa(mu1,eta1);
    Phi_mu2_eta1 = aa(mu2,eta1);

    Phi_mu1_eta1_PI = (Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1';    % Pseudo-inverse
    q_out_eta1 = Phi_mu1_eta1_PI*y;

    y_mu2_est = Phi_mu2_eta1*q_out_eta1;    % Estimated output
end

%% RMSE

mu1 = out_dof;              % Observed nodes
mu2 = 1:dof; mu2(mu1)=[];   % Unobserved nodes

if opt.condimp == 1
    y_acc = y_acc(:,2:end);
    Y_acc = y_acc(:);
    t = t(2:end);
end

% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:(dof-ms)
    if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:) - y_est(mu2(i),:)).^2))); end
    if opt.method == 'ME'; RMSE(i) = sqrt(mean((y_acc(mu2(i),:) - y_mu2_est(i,:)).^2)); end
end
RMSE_tot = mean(RMSE);




end
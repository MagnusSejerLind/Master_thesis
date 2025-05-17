clc,clear
close all
set(0,'defaultTextInterpreter','latex');
rng('default')

%% System properties & options
opt.sysType = "frame";  % ["chain" / "frame"] - Type of system
opt.method = "ME";      % ["TA"/"ME"] - Virtuel sensing method (Toeplitz's/Modal expansion)
opt.out_type = 0;       % [disp=0 / vel=1 / acc=2] - Define output type
opt.error_mod = 1;      % [0/1] - Include error modeling and noise
opt.nonlinear = 1;      % [0/1] - Include nonlinearties in the system
opt.nonlinType = 1;     % [0=constant / 1=varied] - Define type of nonlineaties
opt.numDOF = 4;         % [-int.-] - Number of DOF --ONLY FOR CHAIN SYSTEM
opt.psLoads = 1;        % [1/0] - Apply pseodu loads to convert nonlinear system to linear model
opt.plot = 1;           % [0/1] - plots results
opt.animate = 0;        % [0/1] - Animates the displacements of the structure
opt.aniSave = 0;        % [0/1] - Save animation
opt

in_dof = [1 3];         % Input DOF
% out_dof = [1 2];        % Output DOF
out_dof = [1 2 3 4 5 6 7 8 9 10 12 15 16 18 19 20 22 23 24];        % Output DOF
% frame dof: (x,y,θ)
%% System modeling

[dof,m,k,xi] = systemSetup(opt);
r = numel(in_dof);
ms = numel(out_dof);

% IC
d0 = zeros(dof,1);
v0 = zeros(dof,1);
z0 = [d0;v0];

% Time
N = 100;
dt = 0.01;
t = 0:dt:(N-1)*dt;

% Input (dofs defined earlier)
u_mag = 100;
u = ones(r,N)*u_mag;
u = u.*sin(t*5);
% u = zeros(r,N);
% u(N*0.2:N*0.3) = u_mag;
% u = u.*rand(size(u));
U = u(:);
if opt.error_mod == 1
    [~,~,snr] = modeling_error(0,0);
    U = awgn(U,snr,'measured');
    u = reshape(U,r,N);
end


% Actucal system - no mod. error, no noise
if opt.sysType == "chain"; [M_acc,~,K_acc] = chain(m,m*0,k,dof); end
addBeamError = 0;
if opt.sysType == "frame"; [M_acc,K_acc,dof,snr] = beamStruc(opt,addBeamError); end
[Phi_acc,Lambda_acc] = eig(K_acc,M_acc);    % modal and spectral matrix
[omegaN_acc,i2] = sort(sqrt(diag(Lambda_acc))); % Natural freq.
omegaN_acc = real(omegaN_acc);
Phi_acc = Phi_acc(:,i2);
dd_acc = sqrt(diag(Phi_acc'*M_acc*Phi_acc));
aa_acc = Phi_acc*diag(1./dd_acc);    % Mass-normalized Phi (eigenvec.)
[alpha_acc,beta_acc] = raylieghDamp(omegaN_acc,xi);
C_acc = alpha_acc*M_acc + beta_acc*K_acc;
C_modal_acc = round(Phi_acc'*C_acc*Phi_acc,10);


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
[alpha,beta] = raylieghDamp(omegaN,xi);
C = alpha*M + beta*K;
C_modal = round(Phi'*C*Phi,10);
xi = diag(C_modal) ./ (2*omegaN);

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
    if opt.nonlinType == 0
    % varied / const. mag.
        cf_nl = 0.1;    % coeffcient of nonlinear damping
        kf_nl = 0.1;    % coeffcient of nonlinear stiffness
    else
        cf_nl = rand(1,dof)*0.1;
        kf_nl = rand(1,dof)*0.1;
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
for i = 1:N
    fd_nl(dof+1:end) = cf_nl*z_old(dof+1:end).*abs(z_old(dof+1:end));   % non-linear damping force (velocity dependent)
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);                        % non-linear stiffness force (displacement dependent)

    z_new = Ad*z_old + Bd*u(:,i) - fd_nl - fk_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);
if opt.error_mod == 1
    Y = awgn(Y,snr,'measured');
    y = reshape(Y,ms,N);
end

% Actual system
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));

fd_nl_acc = zeros(size(z_old_acc));
fk_nl_acc = zeros(size(z_old_acc));
for i = 1:N
    fd_nl_acc(dof+1:end) = cf_nl*z_old_acc(dof+1:end).*abs(z_old_acc(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl_acc(dof+1:end)  = kf_nl*(z_old_acc(1:dof).^3);                            % non-linear stiffness force (displacement dependt)

    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i) - fd_nl_acc + fk_nl_acc;
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);

%% Pseudo loads

if opt.psLoads == 1

    % Full input, org. out
    in_dof_FI = 1:1:dof;
    r_FI = numel(in_dof_FI);
    [Ad_FI,Bd_FI,Cd_FI,Dd_FI] = systemMatriciesSS_dis(M,K,C,dof,in_dof_FI,out_dof,opt.out_type,dt);
    [H_FI] = ToeplitzMatrix(N,ms,r_FI,Ad_FI,Bd_FI,Cd_FI,Dd_FI);

    Gamma = pinv(H_FI)*(Y - H*U);
else
    Gamma = zeros(dof*length(t),1);
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

    Y_est = H_ex*U + H_FIFO*Gamma;

    y_est = reshape(Y_est, dof, N);  % decollapse dof columns

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
%% Visualization

% figure()
% tiledlayout('flow')
% for i = 1:dof
%     nexttile
%     hold on
%     plot(t,Y_acc(i:dof:end),'--')
%     plot(t,Y_est(i:dof:end))
%     legend('actual','est.')
%     title(sprintf('DOF: %d', i));
% end



mu1 = out_dof;   % Observed nodes {y = y_ex(mu1,:)}
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes

if opt.plot == 1
    figure()
    tiledlayout('flow')
    if opt.method == 'TA'; sgtitle("Output estimation - Toeplitz's approach",'Interpreter','latex'); end
    if opt.method == "ME"; sgtitle("Output estimation - Modal expansion",'Interpreter','latex'); end

    for i = 1:numel(mu2)
        nexttile
        plot(t,y_acc(mu2(i),:)','k',LineWidth=2)
        hold on

        if opt.method == 'TA'; plot(t,y_est(mu2(i),:),'r--',LineWidth=2); end
        if opt.method == "ME"; plot(t,y_mu2_est(i,:),'--r',LineWidth=2); end

        legend('Actual output','Estimated output')
        title(sprintf('DOF: %d', mu2(i)));
        grid
        xlabel('Time [s]')
        ylabel(sprintf('Output (%d)', opt.out_type));
        xlim([0 N*dt])
    end
end

%% RMSE

% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:(dof-ms)
    if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:) - y_est(mu2(i),:)).^2))); end
    if opt.method == 'ME'; RMSE(i) = sqrt(mean((y_acc(mu2(i),:) - y_mu2_est(i,:)).^2)); end
end
RMSE_tot = mean(RMSE)


%%
% Animate displacements
if opt.animate == 1 && opt.out_type == 0 && opt.sysType == "frame"
    beamStrucDisplacement(y_est,u,in_dof,opt)
end

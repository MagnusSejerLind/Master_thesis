clc,clear,
close all
set(0,'defaultTextInterpreter','latex');
rng('default')

opt.plot = 1;           % [0/1] - plots results
opt.animate = 0;        % [0/1] - Animates the displacements of the structure
opt.aniSave = 0;        % [0/1] - Save animation
%% System properties
opt.sysType = "chain";  % ["chain" / "frame"] - Type of system
opt.method = "ME";      % ["TA"/"ME"] - Virtuel sensing method (Toeplitz's/Modal expansion)
opt.out_type = 0;       % [disp=0 / vel=1 / acc=2] - Define output type
opt.error_mod = 1;      % [0/1] - Include error modeling and noise
opt.nonlinear = 1;      % [0/1] - Include nonlinearties in the system
opt.nonlinType = 1;     % [0=constant / 1=varied] - Define type of nonlineaties
opt.numDOF = 4;         % [-int.-] - Number of DOF --ONLY FOR CHAIN SYSTEM
opt.psLoads = 0;        % [1/0] - Apply pseodu loads to convert nonlinear system to linear model
opt

in_dof = [1 3];         % Input DOF
out_dof = [1 3];        % Output DOF
% out_dof = [1 2 3 4 5 6 7 8 9 10 12 15 16 18 19 20 22 23 24];        % Output DOF
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
N = 500;
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


% Actucal system - no mod. error
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

%%

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
    fd_nl(dof+1:end) = cf_nl*z_old(dof+1:end).*abs(z_old(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);                        % non-linear stiffness force (displacement dependt)

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

%% Output estimation

% Toeplitz's approach
if opt.method == 'TA'

    Psi = H_ex*pinv(H)*Y;

    % decollapse dof columns
    psi = reshape(Psi, dof_ex, N)';

end


% Modal expansion
if opt.method == 'ME'

    mu1 = out_dof;   % Observed nodes {y = y_ex(mu1,:)}
    mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes
    eta1 = 1:numel(mu1);  % Retained modes

    Phi_mu1_eta1 = aa(mu1,eta1);
    Phi_mu2_eta1 = aa(mu2,eta1);

    Phi_mu1_eta1_PI = (Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1';    % Pseudo-inverse
    q_out_eta1 = Phi_mu1_eta1_PI*y;

    y_mu2_est = Phi_mu2_eta1*q_out_eta1;    % Estimated output
    
end



%% Pseudo load transformation

if opt.psLoads == 1 && opt.nonlinear == 1
    % Calculates pseudo loads to convert the nonlinear system to a linear model

    % Full input, unchanged output config., r=n, m=m
    in_dof_con = 1:1:dof;
    r_con=numel(in_dof_con);

    ms_con = numel(out_dof_ex);  % - full output
    [Ad_con,Bd_con,Cd_con,Dd_con] = systemMatriciesSS_dis(M,K,C,dof,in_dof_con,out_dof_ex,opt.out_type,dt);
    [H_con] = ToeplitzMatrix(N,ms_con,r_con,Ad_con,Bd_con,Cd_con,Dd_con);

    % [Ad_con,Bd_con,Cd_con,Dd_con] = systemMatriciesSS_dis(M,K,C,dof,in_dof_con,out_dof,opt.out_type,dt);
    % [H_con] = ToeplitzMatrix(N,ms,r_con,Ad_con,Bd_con,Cd_con,Dd_con);
    % Gamma = pinv(H_con)*(Y - H*U);  % pseudo loads - using Y at output locations.
    % Y_con = H*U + H_con*Gamma;      % LTI model + pseudo loads term

    Gamma = pinv(H_con)*(Psi - H_ex*U);  % pseudo loads - using Y at output locations.
    Y_con = H_ex*U + H_con*Gamma;      % LTI model + pseudo loads term

    Y = Y_con;
    % y = reshape(Y, ms, N);  % decollapse dof columns
    y = reshape(Y, dof, N);  % decollapse dof columns
end



%% Visualization of estimated output

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

        % if opt.method == 'TA'; plot(t,psi(:,mu2(i)),'r--',LineWidth=2); end
        if opt.method == 'TA'; plot(t,y(mu2(i),:),'r--',LineWidth=2); end
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

% {psi = y} - conv.

% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:(dof-ms)
    % if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:)' - psi(:,mu2(i))).^2))); end
    if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:) - y(mu2(i),:)).^2))); end
    if opt.method == 'ME'; RMSE(i) = sqrt(mean((y_acc(mu2(i),:) - y_mu2_est(i,:)).^2)); end
end
RMSE_tot = mean(RMSE)


%%
% Animate displacements
if opt.animate == 1 && opt.out_type == 0 && opt.sysType == "frame"
    beamStrucDisplacement(y,u,in_dof,opt)
end

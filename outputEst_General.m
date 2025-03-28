clc,clear,
close all
set(0,'defaultTextInterpreter','latex');
rng('default')
%% System properties

opt.sysType = "chain";  % ["chain"]
opt.method = "ME";      % ["TA"/"ME"]
opt.error_mod = 0;      % [0/1]
opt.nonlinear = 1;      % [0/1]
opt.nonlinType = 1;     % [0=constant / 1=varied]
opt.out_type = 0;       % [disp=0 / vel=1 / acc=2]
opt.numDOF = 8
in_dof = [1 3];
out_dof = [1 3];
opt.plot = 1;           % [0/1]
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
u_mag = 10;
u = ones(r,N)*u_mag;
% u = u.*sin(t*10);
% u = zeros(r,N);
% u(N*0.2) = u_mag;

% Actucal system
[M_acc,~,K_acc] = chain(m,m*0,k,dof);
[Phi_acc,Lambda_acc] = eig(K_acc,M_acc);    % modal and spectral matrix
[omegaN_acc,i2] = sort(sqrt(diag(Lambda_acc))); % Natural freq.
omegaN_acc = real(omegaN_acc);
Phi_acc=Phi_acc(:,i2);
dd_acc = sqrt(diag(Phi_acc'*M_acc*Phi_acc));
aa_acc = Phi_acc*diag(1./dd_acc);    % Mass-normalized Phi (eigenvec.)
C_modal_acc = diag(2*xi.*omegaN_acc);
C_acc = inv(aa_acc)'*C_modal_acc*inv(aa_acc);

% Base system
if opt.error_mod == 1; [k,m,snr] = modeling_error(k,m); end
[M,~,K] = chain(m,m*0,k,dof);
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);

% Extended system
in_dof_ex = in_dof;
out_dof_ex = (1:1:dof); % r=n
dof_ex = numel(out_dof_ex);
r_ex=numel(in_dof_ex);
ms_ex=numel(out_dof_ex);

% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,opt.out_type,dt);
[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);

% Toeplitz's matricies
[H] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);
[H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

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
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);    % non-linear stiffness force (displacement dependt)

    z_new = Ad*z_old + Bd*u(:,i) - fd_nl - fk_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);
if opt.error_mod == 1
    Y = awgn(Y,snr,'measured');
end

% Actual system
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));

fd_nl_acc = zeros(size(z_old_acc));
fk_nl_acc = zeros(size(z_old_acc));
for i = 1:N
    fd_nl_acc(dof+1:end) = cf_nl*z_old_acc(dof+1:end).*abs(z_old_acc(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl_acc(dof+1:end)  = kf_nl*(z_old_acc(1:dof).^3);    % non-linear stiffness force (displacement dependt)

    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i) - fd_nl_acc + fk_nl_acc;
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);

%% Estimated output

% Teoplitz's approach
if opt.method == 'TA'

    Gamma = H_ex*pinv(H)*Y;

    % inv. dof columns
    gamma = zeros(N,dof_ex);
    for i = 1:dof_ex
        gamma(:,i) = Gamma(i:dof_ex:end);
    end
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
%% Visualization of estimated output

mu1 = out_dof;   % Observed nodes {y = y_ex(mu1,:)}
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes

if opt.plot == 1
    figure()
    tiledlayout('flow')
    if opt.method == 'TA'; sgtitle("Output estimation - Teoplitz's approach",'Interpreter','latex'); end
    if opt.method == "ME"; sgtitle("Output estimation - Modal expansion",'Interpreter','latex'); end

  
    for i = 1:numel(mu2)
        nexttile
        plot(t,y_acc(mu2(i),:)','k',LineWidth=2)
        hold on

        if opt.method == 'TA'; plot(t,gamma(:,mu2(i)),'r--',LineWidth=2); end
        if opt.method == "ME"; plot(t,y_mu2_est(i,:),'--r',LineWidth=2); end

        legend('Actual output','Estimated output')
        title(sprintf('DOF: %d', mu2(i)));
        grid
        xlabel('Time [s]')
        ylabel(sprintf('Output (%d)', opt.out_type));
    end
end

%% Difference

% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:(dof-ms)
    if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:)' - gamma(:,mu2(i))).^2))); end
    if opt.method == 'ME'; RMSE(i) = sqrt(mean((y_acc(mu2(i),:) - y_mu2_est(i,:)).^2)); end
end
RMSE_tot = mean(RMSE)


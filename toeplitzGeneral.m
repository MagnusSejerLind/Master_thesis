clc,clear,close all
set(0,'defaultTextInterpreter','latex');

%% System properties

opt.sysType = "chain";  % ["chain"]
opt.error_mod = 1; % [0/1]
opt.nonlinear = 1;  % [0/1]
out_type = 0;   % [disp=0, vel=1, acc=2]
in_dof = [1 3];
out_dof = [1 3];



%% System modeling

[dof,m,k,xi] = systemSetup(opt.sysType);

dt = 0.01;
r=numel(in_dof);
ms=numel(out_dof);

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

% System
if opt.error_mod == 1; [k,m,snr] = modeling_error(k,m); end
[M,~,K] = chain(m,m*0,k,dof);
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);   % Damping matrix

% Extended system
in_dof_ex = in_dof;
out_dof_ex = [1 2 3 4];
dof_ex = numel(out_dof_ex);
[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,out_type,dt);
r_ex=numel(in_dof_ex);
ms_ex=numel(out_dof_ex);

% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,out_type,dt);
[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof_ex,out_dof_ex,out_type,dt);

% IC
d0=ones(dof,1)*0;
v0=ones(dof,1)*0;
z0=[d0 ; v0];

% Time
N = 500;
t = 0:dt:(N-1)*dt;

% Input (nodes defined earlier)
u_mag = 100;
% u = ones(r,N)*u_mag;
% u = u.*sin(t*10);
u = zeros(r,N);
u(N*0.2) = u_mag;


[H] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);
[H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

%% Compute outputs (system)

if opt.nonlinear == 1
    cf_nl = 0.1;    % coeffcient of nonlinear damping force
    kf_nl = 0.1;    % coeffcient of nonlinear stiffness force
else
    cf_nl = 0; 
    kf_nl = 0; 
end

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

if snr ~= 'none' 
Y = awgn(Y,snr,'measured');
end


%% Estimated output
Gamma = H_ex*pinv(H)*Y;

% inv. dof columns
gamma = zeros(N,dof_ex);
for i = 1:dof_ex
    gamma(:,i) = Gamma(i:dof_ex:end);
end

%% Compute actual output


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

%% Visualization of estimated output

% Estimated dof to vis.
mu1 = out_dof;   % Observed nodes 
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes, rest of nodes

dof_est = 4;    

gamma_est = gamma(:,mu2);


figure()
tiledlayout('flow')
sgtitle('Output estimation - Teoplitz approach') 
for i = 1:numel(mu2)
    nexttile
    plot(t,y_acc(mu2(i),:)','k',LineWidth=2)
    hold on
    plot(t,gamma(:,mu2(i)),'r--',LineWidth=2)
    legend('Actual output','Estimated output')
    title(sprintf('dof no.: %d', mu(i)));
    grid
    xlabel('Time [s]')
    ylabel(sprintf('Output (%d)', out_type));
end


%% Difference


% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:ms
RMSE(i) = sqrt(mean((y_acc(mu2(i),:)' - gamma(:,mu2(i))).^2));
end
RMSE_tot = sum(RMSE)





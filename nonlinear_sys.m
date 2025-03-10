
clc,clear
% close all
set(0,'defaultTextInterpreter','latex');

%% LTI system 

sysType = "chain";
[dof,m,k,xi] = systemSetup(sysType);

out_dof = [1 2];
out_type = 0;   % disp=0, vel=1, acc=2
in_dof = [1 2];
dt = 0.01;
r=numel(in_dof);
ms=numel(out_dof);

% Extended (actual system)
in_dof_ex = in_dof;
out_dof_ex = [1 2 3 4];
dof_ex = numel(out_dof_ex);
r_ex=numel(in_dof_ex);
ms_ex=numel(out_dof_ex);


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

% Model errors
[k,m] = modeling_error(k,m);
[M,~,K] = chain(m,m*0,k,dof);
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);   % Damping matrix

% System matricies S-S
[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof_ex,out_dof_ex,out_type,dt);
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,out_type,dt);

% IC
d0=ones(dof,1)*0;
v0=ones(dof,1)*1;
z0=[d0 ; v0];

% Time
N = 500;
t = 0:dt:(N-1)*dt;

% Input (nodes defined earlier)
u_mag = 100;
u = ones(r,N)*u_mag;
% u = u.*sin(t*10);
% u = zeros(r,N);
% u(N*0.2) = u_mag;
U = u(:);

[H_N_LTI] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd); % (H_N tilde)

%% Compute outputs actual system (Not physical system)
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));
y_acc = zeros(dof,N);
for i = 1:N
    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i);
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);


%% Non-linear system
% (non-extended)

cf_nl = 0.01;    % coeffcient of nonlinear damping force


z_old = z0;
z_new = zeros(size(z_old));


fd_nl = zeros(size(z_old));
for i = 1:N
    fd_nl(dof+1:end) = cf_nl*z_old(dof+1:end).*abs(z_old(dof+1:end));   % non-linear damping force (velocity dependt)
    z_new = Ad*z_old + Bd*u(:,i) - fd_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);

%%

figure()
plot(t,Y(1:r:end),LineWidth=2)
hold on
plot(t,Y_acc(1:dof:end),LineWidth=2)
legend('non-linear','LTI')
grid

%% nonlinear extended system


[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof_ex,out_dof_ex,out_type,dt);

z_old_ex = z0;
z_new_ex = zeros(size(z_old));


fd_nl_ex = zeros(size(z_old));
for i = 1:N
    fd_nl_ex(dof+1:end) = cf_nl*z_old_ex(dof+1:end).*abs(z_old_ex(dof+1:end));   % non-linear damping force (velocity dependt)
    z_new_ex = Ad_ex*z_old_ex + Bd_ex*u(:,i) - fd_nl_ex;
    y_ex(:,i) = Cd_ex*z_old_ex + Dd_ex*u(:,i);
    z_old_ex = z_new_ex;
end
Y_ex = y_ex(:);




clc,clear,close all

set(0,'defaultTextInterpreter','latex');
%% System
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
v0=ones(dof,1)*0;
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


%% Compute outputs actual system

z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));
for i = 1:N
    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i);
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);


%% Modal expansion


mu1 = out_dof;   % Observed nodes 
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes, rest of nodes


eta1 = 1:numel(mu1);  % Retained modes - m=p: determined system

Phi_mu1_eta1 = aa(mu1,eta1);
Phi_mu2_eta1 = aa(mu2,eta1);

Phi_mu1_eta1_MPpi = (Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1';  % Moore-Penrose psedo inversion
acc_mu1 = y_acc(eta1,:);    % Observed acceleration (sim)

q_acc_eta1 = Phi_mu1_eta1_MPpi*acc_mu1;
q_acc_eta1=(Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1'*y_acc(mu1,:);

% modal expansion acceleration estimation
y_mu2_est = Phi_mu2_eta1*q_acc_eta1;

%% Visualization

dof_est = 4;    

figure()
plot(t,y_acc(dof_est,:),'k',LineWidth=2)
hold on
plot(t,y_mu2_est(dof_est-ms,:),'--r',LineWidth=2)
legend('Actual system')
legend('Actual output','Estimated output')
title('Output estimation - Modal expansion')
subtitle(sprintf('dof no.: %d', dof_est));
grid
xlabel('Time [s]')
ylabel(sprintf('Output (%d)', out_type));






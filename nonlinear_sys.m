
clc,clear
close all
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
v0=ones(dof,1)*0;
z0=[d0 ; v0];

% Time
N = 500;
t = 0:dt:(N-1)*dt;

% Input (nodes defined earlier)
u_mag = 100;
% u = ones(r,N)*u_mag;
% u = u.*sin(t*10);
% u = zeros(r,N);
u = ones(r,N);
u(N*0.2) = u_mag;
U = u(:);


%% Compute outputs actual system (extended LTI)
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));
y_acc = zeros(dof,N);
for i = 1:N
    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i);
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);


%% LTI non-extended

cf_nl = 0;    % coeffcient of nonlinear damping force
kf_nl = 0;    % coeffcient of nonlinear stiffness force

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


%% nonlinear extended system

cf_nl = 0.01;    % coeffcient of nonlinear damping force
kf_nl = 0.1;    % coeffcient of nonlinear stiffness force


[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof_ex,out_dof_ex,out_type,dt);

z_old_ex = z0;
z_new_ex = zeros(size(z_old));

fd_nl_ex = zeros(size(z_old));
fk_nl_ex = zeros(size(z_old));
for i = 1:N
    fd_nl_ex(dof+1:end) = cf_nl*z_old_ex(dof+1:end).*abs(z_old_ex(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl_ex(dof+1:end)  = kf_nl*(z_old_ex(1:dof).^3);    % non-linear stiffness force (displacement dependt)

    z_new_ex = Ad_ex*z_old_ex + Bd_ex*u(:,i) - fd_nl_ex + fk_nl_ex;
    y_ex(:,i) = Cd_ex*z_old_ex + Dd_ex*u(:,i);
    z_old_ex = z_new_ex;
end
Y_ex = y_ex(:);


%% Compare LTI to non-linear

figure() 
plot(t,Y_ex(1+1:dof:end),LineWidth=2)
hold on
plot(t,Y_acc(1+1:dof:end),LineWidth=2)
legend('non-linear','LTI')
grid
title('Output: linear/nonlinear')


%% Conversion from physical to LTI

% system models: 
%   Y_ex: nonlinear extended (ms=4,r=2), Y_acc: LTI extended, Y: LTI

    
[H] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);  % (H tilde)
% [H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);  

%% System: m output, n input (H hat): _con {converter}

[dof,m_con,k_con,xi_con] = systemSetup(sysType);

out_dof_con = out_dof;
in_dof_con = 1:1:dof;
r_con=numel(in_dof_con);
ms_con=numel(out_dof_con);

[k_con,m_con] = modeling_error(k_con,m_con);

[M_con,~,K_con] = chain(m_con,m_con*0,k_con,dof);
[Phi_con,Lambda_con] = eig(K_con,M_con);    % modal and spectral matrix
[omegaN_con,i2] = sort(sqrt(diag(Lambda_con))); % Natural freq.
omegaN_con = real(omegaN_con);
Phi_con=Phi_con(:,i2);
dd_con = sqrt(diag(Phi_con'*M_con*Phi_con));
aa_con = Phi_con*diag(1./dd_con);    % Mass-normalized Phi (eigenvec.)
C_modal_con = diag(2*xi_con.*omegaN_con);
C_con = inv(aa_con)'*C_modal_con*inv(aa_con);   % Damping matrix

[Ad_con,Bd_con,Cd_con,Dd_con] = systemMatriciesSS_dis(M_con,K_con,C_con,dof,in_dof_con,out_dof_con,out_type,dt);

[H_con] = TeoplitzMatrix(N,ms_con,r_con,Ad_con,Bd_con,Cd_con,Dd_con);


%% 

% Y_unknown_initial = zeros(ms*N,1);  
% H_con_inv = pinv(H_con);
% 
% L2_Gamma = @(Y_unknown) norm(H_con_inv * (Y_unknown - H * U),2)^2;  % L2 norm of Gamma
% 
% options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton','MaxIterations', 100,'TolX', 1e-10);  % Optimization options
% 
% % Minimize the L2 norm 
% Y_opt = fminunc(L2_Gamma, Y_unknown_initial, options);
% 
% L2_Gamma_opt =  norm(H_con_inv * (Y_opt - H * U),2)^2;  % L2 norm of Gamma
% 
% Theta_opt = Y_opt - H * U;
% Gamma_opt = pinv(H_con) * Theta_opt;
% Y_calc = H*U + H_con*Gamma_opt;
% %% Y=Y_ex

% Y_unknown_initial = zeros(ms*N,1);  
% H_con_inv = pinv(H_con);
% 
% L2_Gamma = @(Y_unknown) norm(H_con_inv * (Y_unknown - H * U),2)^2;  % L2 norm of Gamma
% 
% options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton','MaxIterations', 100,'TolX', 1e-10);  % Optimization options
% 
% % Minimize the L2 norm 
% Y_opt = fminunc(L2_Gamma, Y_unknown_initial, options);
% 
% L2_Gamma_opt =  norm(H_con_inv * (Y_opt - H * U),2)^2;  % L2 norm of Gamma
% 
% Theta_opt = Y_opt - H * U;
% Gamma_opt = pinv(H_con) * Theta_opt;
% Y_calc = H*U + H_con*Gamma_opt;



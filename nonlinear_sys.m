
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
v0=ones(dof,1)*1;
z0=[d0 ; v0];

% Time
N = 500;
t = 0:dt:(N-1)*dt;

% Input (nodes defined earlier)
u_mag = 100;
u = ones(r,N)*u_mag;
u = u.*sin(t*10);
% u = zeros(r,N);
% u = ones(r,N);
% u(N*0.2) = u_mag;
U = u(:);

% [H_N_LTI] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd); % (H_N tilde)

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
z_old = z0;
z_new = zeros(size(z_old));

for i = 1:N
    z_new = Ad*z_old + Bd*u(:,i);
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);

%%

% figure() - no longer correct
% plot(t,Y(1:r:end),LineWidth=2)
% hold on
% plot(t,Y_acc(1:dof:end),LineWidth=2)
% legend('non-linear','LTI')
% grid
% title('Output: linear/nonlinear')

%% nonlinear extended system

cf_nl = 0.1;    % coeffcient of nonlinear damping force

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


%%

% 
% figure() 
% plot(t,Y_ex(1:dof:end),LineWidth=2)
% hold on
% plot(t,Y_acc(1:dof:end),LineWidth=2)
% legend('non-linear','LTI')
% grid
% title('Output: linear/nonlinear')


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


%% For known outout
% Y = H*U+H_ex*Gamma
% Gamma = pinv(H_ex).*Y - pinv(H_ex).*H.*U
% Theta = Y - H*U;
% Gamma = H_ex*Theta;

% 
% Y_conv = Y_acc - Gamma;
% 
% figure()
% hold on
% plot(t,Y_acc(1:dof:end),'k',LineWidth=2)
% plot(t,Y_conv(1:dof:end),'--r',LineWidth=2)
% plot(t,Y(1:r:end))
% legend('LTI model','Converted form phys','non-linear')
% grid
% xlabel('Time [s]')
% ylabel('Output')
% title('Known input')

%% Unknown output


% % L2 norm
% 
% L2 = sqrt(sum(Gamma.^2))
% 
% Theta = Y_unknown - H*U;
% Gamma = H_ex*Theta;


%%

H_con_inv = pinv(H_con);

% Define your objective function (L2 norm)
objective = @(Y_unknown) norm(H_con_inv * (Y_unknown - H * U))^2;
Y_unknown_initial = zeros(ms*N,1);  
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton','MaxIterations', 10,'TolX', 1e-9);  % Optimization options

% Minimize the L2 norm 
Y_opt = fminunc(objective, Y_unknown_initial, options);

% Compute the minimized L2 norm
L2_minimized = sqrt(sum((pinv(H_con) * (Y_opt - H * U)).^2));


%%

Theta_opt = Y_opt - H * U;
Gamma_opt = pinv(H_con) * Theta_opt;

Y_calc = H*U + H_con*Gamma_opt;



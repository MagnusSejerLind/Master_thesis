clc,clear
% close all

set(0,'defaultTextInterpreter','latex');

%% System properties

sysType = "chain";
[dof,m,k,xi] = systemSetup(sysType);

out_dof = [1 3];
out_type = 0;   % disp=0, vel=1, acc=2
in_dof = [1 3];
dt = 0.01;
r=numel(in_dof);
ms=numel(out_dof);

%% System modeling

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
[k,m,snr] = modeling_error(k,m);


[M,~,K] = chain(m,m*0,k,dof);
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);   % Damping matrix

% System matricies
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

%% Compute outputs

y=zeros(ms,N);
z_old=z0;

z_new = zeros(size(z_old));
for i = 1:N
    z_new = Ad*z_old + Bd*u(:,i);
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end

% Output visualization
% plot(t,y(numel(out_dof),:))
% plot(t,y)

%% Teoplitz matrix
[H_N] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);


%% Check consistency

% Rearragne into column vectors 
Y = y(:);
U = u(:);

Y_check = H_N*U;
U_check = pinv(H_N)*Y;
U_consis = U;

% If output is not acceleration, does H_N not have full rank, removes first
% row and last column, and last entry of output vector U (due to no
% information), and Y
if out_type == 0 | out_type == 1
    % H_N_check = H_N(2:end, 1:end-1);
    U_check = U_check(1:end-ms);
    Y_check = Y_check(1:end-ms);
end

% figure()
% plot(U_consis(1:ms:end),'k',LineWidth=2)
% hold on
% plot(U_check(1:ms:end),'r--',LineWidth=2)
% title('Input consistency')
% xlabel('Time')
% ylabel('Input')
% legend('Actual input','Calculated input')
% ylim([min(U_consis)*0.5, max(U_consis)*1.5])
% grid

%% Expanded system

in_dof_ex = in_dof;
out_dof_ex = [1 2 3 4];
dof_ex = numel(out_dof_ex);
[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,out_type,dt);
r_ex=numel(in_dof_ex);
ms_ex=numel(out_dof_ex);

H_N_ex=[];
for jj=1:N
    if jj==1
        H_N_ex((jj-1)*ms_ex+1:jj*ms_ex,1:jj*r_ex)=[Dd_ex H_N_ex];
    else
        Q=Cd_ex*Ad_ex^((jj-1)-1)*Bd_ex;
        H_N_ex((jj-1)*ms_ex+1:jj*ms_ex,1:jj*r_ex)=[Q H_N_ex((jj-2)*ms_ex+1:(jj-1)*ms_ex,:)];
    end
end


%% Estimated expanded output

% output noise, [snr defined in modeling_error function]
if snr ~= 'none' 
Y = awgn(Y,snr,'measured');
end

Gamma = H_N_ex*pinv(H_N)*Y;

% inv. dof columns
gamma = zeros(N,dof_ex);
for i = 1:dof_ex
    gamma(:,i) = Gamma(i:dof_ex:end);
end


%% Expanded system outout  

% y_ex=zeros(dof_ex,N);
% z_old_ex=z0;
% 
% z_new_ex = zeros(size(z_old_ex));
% for i = 1:N
%     z_new_ex = Ad_ex*z_old_ex + Bd_ex*u(:,i);
%     y_ex(:,i) = Cd_ex*z_old_ex + Dd_ex*u(:,i);
%     z_old_ex = z_new_ex;
% end
% Y_ex = y_ex(:);

%% Actual system output

[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof_ex,out_dof_ex,out_type,dt);

z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));
for i = 1:N
    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i);
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);


%% Visualization of estimated output

% Estimated dof
dof_est = 4;    

gamma_est = gamma(:,dof_est);


figure()
% plot(t,y_ex(dof_est,:)','k',LineWidth=2)
plot(t,y_acc(dof_est,:)','k',LineWidth=2)
hold on
plot(t,gamma_est,'r--',LineWidth=2)
legend('Actual output','Estimated output')
title('Output estimation - Teoplitz approch')
subtitle(sprintf('dof no.: %d', dof_est));
grid
xlabel('Time [s]')
ylabel(sprintf('Output (%d)', out_type));


% diff_est = sum(y_acc(dof_est,:)'.^2 - gamma_est.^2);
% % if diff_est < 1e-10; diff_est = 0; end
% text(t(end)*0.7, min(y_acc(dof_est,:))*0.8, ['diff: ', num2str(diff_est)], 'FontSize', 10, 'Color', 'red', 'BackgroundColor', 'white');

%% Difference - Estimated DOFs


mu1 = out_dof;   % Observed nodes 
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes, rest of nodes

% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:ms
RMSE(i) = sqrt(mean((y_acc(mu2(i),:)' - gamma(:,mu2(i))).^2));
end
RMSE_tot = sum(RMSE)



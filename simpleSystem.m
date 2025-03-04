
clc,clear,close all
set(0,'defaultTextInterpreter','latex');
%% System properties
dof = 4;
m = ones(1,dof)*1;
k = ones(1,dof)*300;
xi = (ones(1,dof)*0.1)';
out_dof = [1 2];
out_type = 0;   % disp=0, vel=1, acc=2
in_dof = [1 2];
dt = 0.01;
r=numel(in_dof);
ms=numel(out_dof);

%% Second order modeling

[M,~,K] = chain(m,m*0,k,dof);

[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
Phi=Phi(:,i2);
% Apply mass-normalization
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)

% Damping matrix
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);


%% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,out_type,dt);

%% System

% IC
d0=ones(dof,1)*0;
v0=ones(dof,1)*0;
z0=[d0 ; v0];

% Time
N = 500;
t = 0:dt:(N-1)*dt;

% Input (nodes defined earlier)
u_mag = 1;
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

%% Output visualization

% plot(t,y(numel(out_dof),:))
% plot(t,y)

%% Teoplitz block matrix

% ms: no. of output dof
% r: no. of input dof

H_N=[];
for jj=1:N
    
    if jj==1
        H_N((jj-1)*ms+1:jj*ms,1:jj*r)=[Dd H_N];
    else
        Q=Cd*Ad^((jj-1)-1)*Bd;
        H_N((jj-1)*ms+1:jj*ms,1:jj*r)=[Q H_N((jj-2)*ms+1:(jj-1)*ms,:)];
    end
end


%% Check consistency

% Rearragne into column vectors 
Y = y(:);
U = u(:);


% If output is not acceleration, does H_N not have full rank, removes first
% row and last column, and last entry of output vector U (due to no
% information), and Y
if out_type == 0 | out_type == 1
    H_N_consis = H_N(2:end, 1:end-1);
    U_consis = U(1:end-ms);
    Y_consis = Y(1:end-ms);
end


% Y_check = H_N_consis*U_consis;
% U_check = pinv(H_N_consis)*Y_consis;
Y_check = H_N*U;
U_check = pinv(H_N)*Y;

% y_dof3 = Y(3:dof:end);


if out_type == 0 | out_type == 1
    % H_N_check = H_N(2:end, 1:end-1);
    U_check = U_check(1:end-ms);
    Y_check = Y_check(1:end-ms);
end

figure()
plot(U_consis,'k',LineWidth=2)
hold on
plot(U_check,'r--',LineWidth=2)
title('Input consistency')
xlabel('Time')
ylabel('Input')
legend('Actual input','Calculated input')
ylim([min(U_consis)*0.5, max(U_consis)*1.5])

%% Expanded system


in_dof_ex = in_dof;
out_dof_ex = [1 2 3];
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

% if out_type == 0 | out_type == 1
%     H_N_ex = H_N_ex(2:end, 1:end-1);
% end


%% Estimated expanded output


% Y = awgn(Y,10,'measured');  % add noise to input

Gamma = H_N_ex*pinv(H_N)*Y;


% inv. dof columns
gamma = zeros(N,dof_ex);
for i = 1:dof_ex
    gamma(:,i) = Gamma(i:dof_ex:end);
end


%% Extended system simulation

y_ex=zeros(dof_ex,N);
z_old_ex=z0;

z_new_ex = zeros(size(z_old_ex));
for i = 1:N
    z_new_ex = Ad_ex*z_old_ex + Bd_ex*u(:,i);
    y_ex(:,i) = Cd_ex*z_old_ex + Dd_ex*u(:,i);
    z_old_ex = z_new_ex;
end
Y_ex = y_ex(:);


%% Visualization of estimated output

% Estimated dof
dof_est = 3;    

gamma_est = gamma(:,dof_est);
y_accEst = y_ex(dof_est,:)';


figure()
plot(y_accEst,'k',LineWidth=2)
hold on
plot(gamma_est,'r--',LineWidth=2)
legend('Actual output','Estimated output')
title('Output estimation')
subtitle(sprintf('dof no.: %d', dof_est));
grid
xlabel('Time')
ylabel(sprintf('Output (%d)', out_type));


diff_est = sum(y_accEst.^2 - gamma_est.^2);
if diff_est < 1e-10; diff_est = 0; end
text(N*0.7, min(y_accEst)*0.8, ['diff: ', num2str(diff_est)], 'FontSize', 10, 'Color', 'red', 'BackgroundColor', 'white');


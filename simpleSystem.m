
clc,clear,close all
%% System properties
dof = 4;
m = ones(1,dof)*1;
k = ones(1,dof)*300;
xi = (ones(1,dof)*0.1)';
out_dof = [1 4];
out_type = 0;   % disp=0, vel=1, acc=2
in_dof = [3 4];
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

%% Compute outputs

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
u = u.*sin(t);
% u = zeros(r,N);
% u(N*0.1) = 1;


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



%% Teoplitz block matrix

% ms: # output dof
% r: # input dof

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


Y_check = H_N*U;
U_check = pinv(H_N)*Y;
% y_check = H_N*u';
% u_check = pinv(H_N)*y';


% y_dof3 = Y(3:dof:end);

%% Expanded system

in_dof_ex = [1 2 3 4];
out_dof_ex = [1 2 3 4];
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


%%
% H_N_ex.*pinv(H_N)














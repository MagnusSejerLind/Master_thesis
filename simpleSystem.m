% Express the state space system matricies for a simple chain system

clc,clear,close all
%% System properties
dof = 4;
m = ones(1,dof)*2;
k = ones(1,dof)*300;
xi = (ones(1,dof)*0.1)';
out_dof = [4];
out_type = 0;   % disp=0, vel=1, acc=2
in_dof = [];
dt = 0.01;

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
Ac=[zeros(dof) eye(dof) ; -M\K -M\C];
r=numel(in_dof);
B2=zeros(dof,r);
for jj=1:r
    B2(in_dof(jj),jj)=1;
end
Bc=[zeros(dof,r) ; M\B2];
ms=numel(out_dof);
if out_type==0
    Cc=[eye(dof) zeros(dof)];
    Cc=Cc(out_dof,:);
    Dc=zeros(ms,r);
elseif out_type==1
    Cc=[zeros(dof) eye(dof)];
    Cc=Cc(out_dof,:);
    Dc=zeros(ms,r);
elseif out_type==2
    Cc=Ac(dof+1:end,:); 
    Cc=Cc(out_dof,:);
    Cacc=eye(dof);
    Cacc=Cacc(out_dof,:);
    Dc=Cacc*M^-1*B2;
else
    disp('Something is wrong in the output type selection.')
    return
end

% Convert to discrete
Ad=expm(Ac*dt);
Bd=Ac^-1*(Ad-eye(2*dof))*Bc;
Cd=Cc;
Dd=Dc;



%% Compute outputs

% IC
d0=ones(dof,1)*1;
v0=ones(dof,1)*0;
z0=[d0 ; v0];

% Time
N = 5000;
t = 0:dt:(N-1)*dt;

u_mag = 0;
u=ones(r,N)*u_mag;

y=zeros(ms,N);
z_old=z0;

z_new = zeros(size(z_old));
for i = 1:N
    z_new = Ad*z_old + Bd*u(:,i);
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end

%%

% output dof1
plot(t,y(1,:))



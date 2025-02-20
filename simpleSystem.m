% Express the state space system matricies for a simple chain system


%% System properties
dof = 4;
m = ones(1,dof)*2;
k = ones(1,dof)*3;
xi = ones(1,4)*0.1;
out_dof = [2 4];
out_type = 0;   % disp=0, vel=1, acc=2
% in_dof = 

%% Second order modeling

[M,~,K] = chain(m,m*0,k,dof);

[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
omegaN = sqrt(diag(Lambda)); % Natural freq.
% accending order
[omegaN,i2]=sort(omegaN);
Phi=Phi(:,i2);
% Apply mass-normalization
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)

% Damping matrix
C_modal = 2*omegaN*xi;
C = inv(aa)'*C_modal*aa;


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










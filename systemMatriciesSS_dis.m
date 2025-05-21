function [Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,out_type,dt)
% Expresses the discrete state space system matricies

% out_type: 0: disp, 1: vel, 2: acc
Ac = [zeros(dof) eye(dof) ; -M\K -M\C];
r = numel(in_dof);

B2 = zeros(dof,r);
for jj = 1:r
    B2(in_dof(jj),jj) = 1;  % input dist. vec.
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

end

% Convert to discrete
Ad = expm(Ac*dt);
Bd = Ac^-1*(Ad-eye(2*dof))*Bc;
Cd = Cc;
Dd = Dc;
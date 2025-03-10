function [H_N] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd)

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
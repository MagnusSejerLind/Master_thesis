function [M,C,K] = chain(m,c,k,dof)
% Expresses system matricies of a chain system


if dof~=1
    for j=1:dof-1
        kd(j) = k(j)+k(j+1);
        cd(j) = c(j)+c(j+1);
        kf(j) = -k(j+1);
        cf(j) = -c(j+1);
    end
    kd(dof) = k(dof);
    cd(dof) = c(dof);
    K = diag(kf,1)+diag(kf,-1)+diag(kd);
    C = (diag(cf,1)+diag(cf,-1)+diag(cd));
    M = diag(m);
else
    K = k;
    M = m;
    C = c;
end
end

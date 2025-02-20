function [M,C,K]=chain(m1,c1,k1,dof)
% Expresses system matricies of a chain system


    if dof~=1
        for j=1:dof-1
           kd(j)=k1(j)+k1(j+1);
           cd(j)=c1(j)+c1(j+1);
           kf(j)=-k1(j+1);
           cf(j)=-c1(j+1);
        end
        kd(dof)=k1(dof);
        cd(dof)=c1(dof);
        K= diag(kf,1)+diag(kf,-1)+diag(kd);
        C=(diag(cf,1)+diag(cf,-1)+diag(cd));
        M=diag(m1);
    else
        K=k1;
        M=m1;
        C=c1;
    end
end
function [ddata,vdata,adata]=newmark_alg(K,C,M,F,alpha,beta,dt,N,d0,v0,ndof)
    % Helping parameters %
    a1=1/(alpha*dt^2);
    a2=1/(alpha*dt);
    a3=1/(2*alpha);
    a4=1/(alpha);
    % Initial accelerations
    a0=M\(F(:,1)-K*d0-C*v0);
    % Effective stiffness matrix %
    keff=a1*M+beta*a2*C+K;
    % Setting up response%
    ddata=zeros(ndof,N); ddata(:,1)=d0;
    vdata=zeros(ndof,N); vdata(:,1)=v0;
    adata=zeros(ndof,N); adata(:,1)=a0;
    % Loop over time %
    for i=1:N-1
        % Effective acceleration term %
        aeff=M*(a1*d0+a2*v0+(a3-1)*a0);
        % Effective velocity term %
        veff=C*(beta*a2*d0+(beta*a4-1)*v0+dt*(beta*a3-1)*a0);    
        % Displacement %
        d=keff\(F(:,i+1)+aeff+veff);
        % Velocity %
        v=beta*a2*(d-d0)-(beta*a4-1)*v0-dt*(beta*a3-1)*a0;
        % Acceleration %
        a=a1*(d-d0-dt*v0)-(a3-1)*a0;    
        % Update of kinematics %
        d0=d;
        v0=v;
        a0=a;      
        % Store responses %
        ddata(:,i+1)=d0;
        vdata(:,i+1)=v0;
        adata(:,i+1)=a0; 
    end
end
function [M, K, dof, snr] = beamStruc(opt,addBeamError)


% Define structure
[file_i, xy, nnod, sizee, idb, dof, incid, l, gamma, m, EA, EJ, T, posiz, nbeam, pr] = loadstructure;


% Draw structure
if opt.plot == 1
    % dis_stru(posiz, l, gamma, xy, pr, idb, dof);
    % ylim([-0.25 2.5])
    % xlim([-0.5 3.5])
end

% Assemble mass and stiffness matrices
[M,K,snr] = assemble(incid, l, m, EA, EJ, T, gamma, idb, opt, addBeamError);  % Obs. modeling error implemented within

% Matrix partitioning - free/constrained
MFF = M(1:dof, 1:dof);
MCF = M(dof+1:end, 1:dof);
MFC = M(1:dof, dof+1:end);
MCC = M(dof+1:end, dof+1:end);

KFF = K(1:dof, 1:dof);
KCF = K(dof+1:end, 1:dof);
KFC = K(1:dof, dof+1:end);
KCC = K(dof+1:end, dof+1:end);

M_full = M;
K_full = K;

% Include only free nodes
M = MFF;
K = KFF;





function [M,K,snr] = assemble(incidenze,l,m,EA,EJ,T,gamma,idb,opt,addBeamError)
% Assembles the mass and stiffness matrix of the beam structure system

[n_el,~]=size(incidenze);
n_dof=max(max(idb));

% Assembling matrices M and K
M=zeros(n_dof,n_dof);
K=zeros(n_dof,n_dof);
for k=1:n_el
    [mG,kG,snr] = elementMat(l(k),m(k),EA(k),EJ(k),T(k),gamma(k),opt,addBeamError);
    for iri=1:6
        for ico=1:6
            i1=incidenze(k,iri);
            i2=incidenze(k,ico);
            M(i1,i2)=M(i1,i2)+mG(iri,ico);
            K(i1,i2)=K(i1,i2)+kG(iri,ico);
        end
    end
end

end

function [mG,kG,snr] = elementMat(l,m,EA,EJ,T,gamma,opt,addBeamError)
% Builds global elemental mass and stiffness matrices


if opt.error_mod == 1 && addBeamError == 1
    [EJ,m,snr] = modeling_error(EJ,m);
else
    snr = "none";
end

% local mass matrix
mL = m*l*[1./3.  0.          0.          1./6.  0.          0.
    0.     13./35.     11.*l/210.  0.     9./70.      -13*l/420.
    0.     11.*l/210.  l^2/105.    0.     13*l/420.   -l^2/140.
    1./6.  0.          0.          1./3.  0.          0.
    0.     9./70.      13*l/420.   0.     13./35.     -11.*l/210.
    0.     -13*l/420.  -l^2/140.   0.     -11.*l/210. l^2/105.   ] ;

% contribution due to the axial deformation
kL_ax = EA/l* [ 1 0 0 -1 0 0
    0 0 0  0 0 0
    0 0 0  0 0 0
    -1 0 0  1 0 0
    0 0 0  0 0 0
    0 0 0  0 0 0 ] ;

% contribution due to bending deformation
kL_fl = EJ * [ 0.    0.       0.      0.    0.       0.
    0.  +12./l^3+T./l./EJ  6./l^2   0.  -12./l^3-T./l./EJ  6./l^2
    0.   6./l^2  +4./l     0.   -6./l^2  +2./l
    0.    0.       0.      0.    0.       0.
    0.  -12./l^3-T./l./EJ  -6./l^2  0.  +12./l^3+T./l./EJ  -6./l^2
    0.   6./l^2  +2./l    0.   -6./l^2  +4./l    ] ;

kL = kL_ax+kL_fl;  % Local stiffness matrix 


% rotation matrix 3x3 (single node DoF)
lambda = [ cos(gamma) sin(gamma) 0.
    -sin(gamma) cos(gamma) 0.
    0.         0.      1.] ;

% rotation matrix 6x6 (element DoF)
Lambda = [ lambda     zeros(3)
    zeros(3)   lambda      ] ;

% Global mass and stiffness matricies
mG = Lambda' * mL * Lambda ;
kG = Lambda' * kL * Lambda ;

end




function dis_stru(posiz,l,gamma,xy,~,idb,dof)
% Plot the structure

xmax = max(xy(:,1));
xmin = min(xy(:,1));
ymax = max(xy(:,2));
ymin = min(xy(:,2));

dx = (xmax - xmin)/100;
dy = (ymax - ymin)/100;
d = sqrt(dx^2 + dy^2);

figure();
hold on;
% elements
for i=1:length(posiz)
    xin=posiz(i,1);
    yin=posiz(i,2);
    xfi=posiz(i,1)+l(i)*cos(gamma(i));
    yfi=posiz(i,2)+l(i)*sin(gamma(i));
    % colore = colori(i,:);
    plot([xin xfi],[yin yfi],'linewidth',2,'color','k');
    %plot([xin xfi],[yin yfi],'b','linewidth',2);
end
grid on; box on;

% nodal positions
plot(xy(:,1),xy(:,2),'r.','markersize',20);

triangolo_h = [ 0 0; -sqrt(3)/2 .5; -sqrt(3)/2 -.5; 0 0]*d*2;
triangolo_v = [ 0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2;
triangolo_r = [0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2 * [sqrt(2)/2 -sqrt(2/2); -sqrt(2)/2 -sqrt(2)/2];

hold on
for ii = 1:size(xy,1)
    %rectangle('Position',[xy(ii,1)-d/2 xy(ii,2)-d/2 d d],'curvature',1,'edgecolor','r','linewidth',3);
    text(xy(ii,1) + d, xy(ii,2) + d,num2str(ii));
    if (idb(ii,1) > dof)
        fill(xy(ii,1) + triangolo_h(:,1),xy(ii,2) + triangolo_h(:,2),'b');
    end
    if (idb(ii,2) > dof)
        fill(xy(ii,1) + triangolo_v(:,1),xy(ii,2) + triangolo_v(:,2),'b');
    end
    if (idb(ii,3) > dof)
        fill(xy(ii,1) + triangolo_r(:,1),xy(ii,2) + triangolo_r(:,2),'b');
    end
end

axis equal
title('Beam Frame Structure')
end


end
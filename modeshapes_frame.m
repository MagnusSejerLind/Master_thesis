clc,clear,close all

set(0,'defaultTextInterpreter','latex');

%% Load input file and assemble structure
[file_i, xy, nnod, sizee, idb, ndof, incid, l, gamma, m, EA, EJ, T, posiz, nbeam, pr] = loadstructure;
    
% %% Draw structure
% dis_stru(posiz, l, gamma, xy, pr, idb, ndof);
% legend('Principal beam element','','','Reinforcement beam element')
% ylim([-0.25 3.5])

%% Check freqeuncy range: max element length

% Examined structure frequency range
freq_range = [0 20];
safetyFactor = 2;

% max length
L_max = sqrt((pi^2)/(safetyFactor*freq_range(2))*min(sqrt(EJ)./m));
if max(l) > L_max
    disp('!!Element length too large!!')
end


%% Assemble mass and stiffness matrices
[M, K] = assem(incid, l, m, EA, EJ, T, gamma, idb);

%% Matrix partitioning

MFF = M(1:ndof, 1:ndof);
MCF = M(ndof+1:end, 1:ndof);
MFC = M(1:ndof, ndof+1:end);
MCC = M(ndof+1:end, ndof+1:end);

KFF = K(1:ndof, 1:ndof);
KCF = K(ndof+1:end, 1:ndof);
KFC = K(1:ndof, ndof+1:end);
KCC = K(ndof+1:end, ndof+1:end);

%% Natural freq. and mode shapes determined

[modes, omega2] = eig(MFF\KFF);
omega = sqrt(diag(omega2));

% Sort frequencies in ascending order
[omega, i_omega] = sort(omega);
freq0 = omega/2/pi;

% Sort mode shapes in ascending order
modes = modes(:, i_omega);

%%
% freq_range = freq_range*1.1;
% Number of natural frequencies within the examined range
freq_in_range = sum(freq0 >= freq_range(1) & freq0 <= freq_range(2));

%% Mode shape plot

scale_factor = 1;

for i = 1:freq_in_range
    mode = modes(:,i);

    figure()
    diseg2(mode,scale_factor,incid,l,gamma,posiz,idb,xy)   
    title(sprintf('Mode %g: f=%.3f Hz, scale factor: %.1f', i,freq0(i),scale_factor))
    legend('Undeformed','Deformed','Interpreter','latex')
    xlabel('x [m]')
    ylabel('y [m]',Rotation=0)
    xlim([-1 4])
    ylim([0 2.5])
end

% 
% % Additional modes
% figure()
% tiledlayout(1,freq_in_range)
% for i = 1:freq_in_range
% 
%     mode = modes(:,i);
% 
%     nexttile
% 
%     diseg2(mode,scale_factor,incid,l,gamma,posiz,idb,xy)   
%     title(sprintf('Mode %g: f=%.3f Hz, scale factor: %.1f', i,freq0(i),scale_factor))
%     legend('Undeformed','Deformed','Interpreter','latex')
%     xlabel('x [m]')
%     ylabel('y [m]',Rotation=0)
%     xlim([-1 5])
%     ylim([0 4.5])
% end



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


function [M,K] = assem(incidenze,l,m,EA,EJ,T,gamma,idb)

% Checking consistency input data
[n_el,nc]=size(incidenze);
if nc ~= 6
    disp('Error: number of columns of incidence matrix different from 6')
    return
end
if length(l) ~= n_el
    sprintf('Error: number of elements in l different from n')
    return
end
if length(m) ~= n_el    
    sprintf('Error: number of elements in m different from number of elements')
    return
end
if length(EA) ~= n_el
    sprintf('Error: number of elements in EA different number of elements')
    return
end
if length(EJ) ~= n_el
    sprintf('Error: number of elements in EJ differenc number of elements')
    return
end
if length(T) ~= n_el
    sprintf('Error: number of elements in T different number of elements')
    return
end
if length(gamma) ~= n_el
    sprintf('Error: number of elements in alpha different number of elements')
    return
end

if min(min(incidenze)) ~= 1    
    sprintf('Error: dof sequence does not start from 1')
    return
end

% Limiting the total dof of the structure
n_dof=max(max(idb));
% if n_dof > 100
%     sprintf('Errore: dof > 100, not allowed!')
%     return
% end

% Assembling matrices M and K
M=zeros(n_dof,n_dof);
K=zeros(n_dof,n_dof);
for k=1:n_el
    [mG,kG] = el_tra(l(k),m(k),EA(k),EJ(k),T(k),gamma(k));
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

function [mG,kG] = el_tra (l,m,EA,EJ,T,gamma)
% mass matrix in the local reference frame
mL = m*l*[1./3.  0.          0.          1./6.  0.          0.
          0.     13./35.     11.*l/210.  0.     9./70.      -13*l/420.
          0.     11.*l/210.  l^2/105.    0.     13*l/420.   -l^2/140.
          1./6.  0.          0.          1./3.  0.          0.
          0.     9./70.      13*l/420.   0.     13./35.     -11.*l/210.
          0.     -13*l/420.  -l^2/140.   0.     -11.*l/210. l^2/105.   ] ;

% stiffness matrix in the local reference frame
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
           
kL = kL_ax+kL_fl ;

% matrix transformation from local to global reference frame
% rotation matrix 3x3 (single node DoF)
lambda = [ cos(gamma) sin(gamma) 0. 
          -sin(gamma) cos(gamma) 0.
              0.         0.      1.] ;

% rotation matrix 6x6 (element DoF)
Lambda = [ lambda     zeros(3)
           zeros(3)   lambda      ] ;

mG = Lambda' * mL * Lambda ;
kG = Lambda' * kL * Lambda ;

end

function diseg2(mode,scale_factor,incidenze,l,gamma,posiz,idb,xy)

% Checking consistency input data
[n_el,nc]=size(incidenze);

if length(posiz) ~= n_el
    sprintf('Error: number of nodes in posit matrix differs from number of elements')
    return
end

n_gdl=length(mode);


hold on
% looping on the finite elements
for k=1:n_el
% building the nodal displacements vector of each element in the global
% reference frame
    for iri=1:6
        if incidenze(k,iri) <= n_gdl
        xkG(iri,1)=mode(incidenze(k,iri));
        else
        xkG(iri,1)=0.;
        end
    end
% Applying scale factor
    xkG=scale_factor*xkG;
% Global to Local reference frame rotation
    lambda = [ cos(gamma(k)) sin(gamma(k)) 0. 
              -sin(gamma(k)) cos(gamma(k)) 0.
	                0.      0.     1. ] ;
    Lambda = [ lambda     zeros(3)
              zeros(3)   lambda      ] ;
    xkL=Lambda*xkG;

% Computing the axial (u) and transversal (w) displacements by means of
% shape functions
    csi=l(k)*[0:0.05:1];
    fu=zeros(6,length(csi));
    fu(1,:)=1-csi/l(k);
    fu(4,:)=csi/l(k);
    u=(fu'*xkL)';
    fw=zeros(6,length(csi));
    fw(2,:)=2*(csi/l(k)).^3-3*(csi/l(k)).^2+1;
    fw(3,:)=l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k));
    fw(5,:)=-2*(csi/l(k)).^3+3*(csi/l(k)).^2;
    fw(6,:)=l(k)*((csi/l(k)).^3-(csi/l(k)).^2);
    w=(fw'*xkL)';
% Local to global transformation of the element's deformation
    xyG=lambda(1:2,1:2)'*[u+csi;w];
    undef=lambda(1:2,1:2)'*[csi;zeros(1,length(csi))];
 % Plotting deformed and undeformed element's shape
    plot(undef(1,:)+posiz(k,1),undef(2,:)+posiz(k,2),'b--')
    plot(xyG(1,:)+posiz(k,1),xyG(2,:)+posiz(k,2),'b')
end

% Looping the nodes
n_nodi=size(idb,1);
xkG=[];
for k=1:n_nodi
    for ixy=1:2
        if idb(k,ixy) <= n_gdl
        xkG(k,ixy)=mode(idb(k,ixy));
        else
        xkG(k,ixy)=0.;
        end
    end
end
xkG=scale_factor*xkG;
xyG=xkG+xy;
plot(xy(:,1),xy(:,2),'b.')
plot(xyG(:,1),xyG(:,2),'bo')

grid on
box on
axis equal
end
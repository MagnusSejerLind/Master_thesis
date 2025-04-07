clc,clear,
% close all
set(0,'defaultTextInterpreter','latex');
rng('default')
opt.plot = 1;           % [0/1]
%% System properties

opt.sysType = "beamStruc";  % ["chain"] - Type of system
opt.method = "TA";      % ["TA"/"ME"] - Virtuel sensing method (Toeplitz's/Modal expansion)
opt.out_type = 2;       % [disp=0 / vel=1 / acc=2] - Define output type
opt.error_mod = 1;      % [0/1] - Include error modeling and noise
opt.nonlinear = 1;      % [0/1] - Include nonlinearties in the system
opt.nonlinType = 1;     % [0=constant / 1=varied] - Define type of nonlineaties


% opt.numDOF = 8;          % Number of DOF

in_dof = [1 3];         % Input DOF
out_dof = [1 2 3 4 5 6 7 8 9 10 12 15 16 18 19 20 22 23 24];        % Output DOF
opt

%% Define structure
[file_i, xy, nnod, sizee, idb, dof, incid, l, gamma, m, EA, EJ, T, posiz, nbeam, pr] = loadstructure;

% Draw structure
if opt.plot == 1
    % dis_stru(posiz, l, gamma, xy, pr, idb, dof);
    % ylim([-0.25 2.5])
    % xlim([-0.5 3.5])
end

% Assemble mass and stiffness matrices
[M,K,snr] = assem(incid, l, m, EA, EJ, T, gamma, idb, opt);

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


%% System modeling

% [dof,m,k,xi] = systemSetup(opt);
xi = (ones(1,dof)*0.01)';

r = numel(in_dof);
ms = numel(out_dof);

% IC
d0 = zeros(dof,1);
v0 = zeros(dof,1);
z0 = [d0;v0];

% Time
N = 500;
dt = 0.01;
t = 0:dt:(N-1)*dt;

% Input (dofs defined earlier)
u_mag = 10;
% u = ones(r,N)*u_mag;
% u = u.*sin(t*5);
u = zeros(r,N);
u(N*0.2) = u_mag;


% Actucal system
M_acc = M;
K_acc = K;
[Phi_acc,Lambda_acc] = eig(K_acc,M_acc);    % modal and spectral matrix
[omegaN_acc,i2] = sort(sqrt(diag(Lambda_acc))); % Natural freq.
omegaN_acc = real(omegaN_acc);
Phi_acc=Phi_acc(:,i2);
dd_acc = sqrt(diag(Phi_acc'*M_acc*Phi_acc));
aa_acc = Phi_acc*diag(1./dd_acc);    % Mass-normalized Phi (eigenvec.)
C_modal_acc = diag(2*xi.*omegaN_acc);
C_acc = inv(aa_acc)'*C_modal_acc*inv(aa_acc);


% Base system
% if opt.error_mod == 1; [k,m,snr] = modeling_error(k,m); end       % Not correct imp.
% M,K already defined
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);


% Extended system
in_dof_ex = in_dof;
out_dof_ex = (1:1:dof); 
% % out_dof_ex = out_dof;
% % in_dof_ex = (1:1:dof); 
dof_ex = numel(out_dof_ex);
r_ex = numel(in_dof_ex);
ms_ex = numel(out_dof_ex);

% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,opt.out_type,dt);
[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof,out_dof_ex,opt.out_type,dt);

% Toeplitz's matricies
[H] = TeoplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);
[H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);


% Nonlinearities
if opt.nonlinear == 1
    if opt.nonlinType == 0
        cf_nl = 0.1;    % coeffcient of nonlinear damping 
        kf_nl = 0.1;    % coeffcient of nonlinear stiffness 
    else
    cf_nl = rand(1,dof)*0.1;    
    kf_nl = rand(1,dof)*0.1;    
    end
else
    cf_nl = 0;
    kf_nl = 0;
end

%% Compute outputs

% Base system
z_old = z0;
z_new = zeros(size(z_old));
fd_nl = zeros(size(z_old));
fk_nl = zeros(size(z_old));
for i = 1:N
    fd_nl(dof+1:end) = cf_nl*z_old(dof+1:end).*abs(z_old(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);    % non-linear stiffness force (displacement dependt)

    z_new = Ad*z_old + Bd*u(:,i) - fd_nl - fk_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);
if opt.error_mod == 1

    Y = awgn(Y,snr,'measured');
end

% Actual system
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));

fd_nl_acc = zeros(size(z_old_acc));
fk_nl_acc = zeros(size(z_old_acc));
for i = 1:N
    fd_nl_acc(dof+1:end) = cf_nl*z_old_acc(dof+1:end).*abs(z_old_acc(dof+1:end));   % non-linear damping force (velocity dependt)
    fk_nl_acc(dof+1:end)  = kf_nl*(z_old_acc(1:dof).^3);    % non-linear stiffness force (displacement dependt)

    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i) - fd_nl_acc + fk_nl_acc;
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);




%% Estimated output

% Teoplitz's approach
if opt.method == 'TA'

    Gamma = H_ex*pinv(H)*Y;

    % inv. dof columns
    gamma = zeros(N,dof_ex);
    for i = 1:dof_ex
        gamma(:,i) = Gamma(i:dof_ex:end);
    end
end


% Modal expansion
if opt.method == 'ME'

    mu1 = out_dof;   % Observed nodes {y = y_ex(mu1,:)}
    mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes   
    eta1 = 1:numel(mu1);  % Retained modes

    Phi_mu1_eta1 = aa(mu1,eta1);
    Phi_mu2_eta1 = aa(mu2,eta1);

    Phi_mu1_eta1_PI = (Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1';    % Pseudo-inverse
    q_out_eta1 = Phi_mu1_eta1_PI*y;

    y_mu2_est = Phi_mu2_eta1*q_out_eta1;    % Estimated output
end



%% Visualization of estimated output

mu1 = out_dof;   % Observed nodes {y = y_ex(mu1,:)}
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes

if opt.plot == 1
    figure()
    tiledlayout('flow')
    if opt.method == 'TA'; sgtitle("Output estimation - Teoplitz's approach",'Interpreter','latex'); end
    if opt.method == "ME"; sgtitle("Output estimation - Modal expansion",'Interpreter','latex'); end

  
    for i = 1:numel(mu2)
        nexttile
        plot(t,y_acc(mu2(i),:)','k',LineWidth=2)
        hold on

        if opt.method == 'TA'; plot(t,gamma(:,mu2(i)),'r--',LineWidth=2); end
        if opt.method == "ME"; plot(t,y_mu2_est(i,:),'--r',LineWidth=2); end

        legend('Actual output','Estimated output')
        title(sprintf('DOF: %d', mu2(i)));
        grid
        xlabel('Time [s]')
        ylabel(sprintf('Output (%d)', opt.out_type));
    end
end


%% Difference

% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:(dof-ms)
    if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:)' - gamma(:,mu2(i))).^2))); end
    if opt.method == 'ME'; RMSE(i) = sqrt(mean((y_acc(mu2(i),:) - y_mu2_est(i,:)).^2)); end
end
RMSE_tot = mean(RMSE)














%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [M,K,snr] = assem(incidenze,l,m,EA,EJ,T,gamma,idb,opt)

[n_el,~]=size(incidenze);

n_dof=max(max(idb));

% Assembling matrices M and K
M=zeros(n_dof,n_dof);
K=zeros(n_dof,n_dof);
for k=1:n_el
    [mG,kG,snr] = el_tra(l(k),m(k),EA(k),EJ(k),T(k),gamma(k),opt);
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

function [mG,kG,snr] = el_tra(l,m,EA,EJ,T,gamma,opt)


%%%%% %%% 
if opt.error_mod == 1
    [EJ,m,snr] = modeling_error(EJ,m);
end
%%%%

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



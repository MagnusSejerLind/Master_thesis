clc,clear,close all

%% Set up system

% Pars
n = 9;  % number of DOFs
L_el = 0.15; % [m] 
L = L_el*(n-1);
A = 0.0001;  % [m^2]
I = 8.3333*1E-10;   % [m^4]
xi_E = 100*1E9;    % [Pa]
rho = 8000; % [kg/m^3]

coor = [(0:L_el:L)'  zeros(n,1)];

%       Node       Node   I       E    A     rho     Element
no = [  1          2      I       xi_E    A     rho     1
        2          3      I       xi_E    A     rho     2
        3          4      I       xi_E    A     rho     3
        4          5      I       xi_E    A     rho     4
        5          6      I       xi_E    A     rho     5
        6          7      I       xi_E    A     rho     6
        7          8      I       xi_E    A     rho     7
        8          9      I       xi_E    A     rho     8];

fixeddof = 1:3;   % x,y,rot at first node

% Express mass and stiffness matricies
[K,M,freedof]=beam2D(coor,no,fixeddof);
dof = numel(freedof);


% Determine eigenfreq. and eigenvec. (mode shapes)
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
omegaN = sqrt(diag(Lambda)); % Natural freq.
% accending order
[omegaN,i2]=sort(omegaN);
Phi=Phi(:,i2);
% Apply mass-normalization
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)

Omega = (omegaN(1)+omegaN(2))/2;


%% Compute damping matrix
xi = [0.02 0.05]';    % Known damping ratio, mode 1,3

% Rayleigh damping
G = [1/(2*omegaN(1)), omegaN(1)/2 ;
    1/(2*omegaN(3)), omegaN(3)/2];
coeff_rayleigh = G^-1*xi;
alpha_rayleigh = coeff_rayleigh(1);
beta_rayleigh = coeff_rayleigh(2);

% Damping matrix
C = alpha_rayleigh*M + beta_rayleigh*K;
% damping ratios
xi_E=diag(aa'*C*aa)./(2*omegaN);


%% Simulate transverse nodal acceleration
% Using time integration (Newmark)

% % Time
N = 2000;   % Signal length
dt = 0.001;
t = 0:dt:(N-1)*dt;

% Integration parameters
alpha = 1/4;
beta = 1/2;

% Initial conditions
d0=zeros(dof,1);
v0=zeros(dof,1);

% Input
F = zeros(dof,N);
F(5,:) = sin(Omega*t);

% Time integration, solving for accelerations  
[~,~,adata]=newmark_alg(K,C,M,F,alpha,beta,dt,N,d0,v0,dof);     % adata: nodal acceleration over time 

%% Modal expansion

mu1 = [5,17,23];   % Observed nodes (transverse direction)
mu2 = 1:dof; mu2(mu1)=[];  % Unobserved nodes, rest of nodes
% node: dof only contains free nodes

eta1 = 1:numel(mu1);  % Retained modes - m=p: determined system

Phi_mu1_eta1 = aa(mu1,eta1);
Phi_mu2_eta1 = aa(mu2,eta1);

Phi_mu1_eta1_MPpi = (Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1';  % Moore-Penrose psedo inversion
acc_mu1 = adata(eta1,:);    % Observed acceleration (from Newmark time integration)

% q_acc_eta1 = Phi_mu1_eta1_MPpi*acc_mu1;
q_acc_eta1=(Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1'*adata(mu1,:);

% modal expansion acceleration estimation
adata_mu2_est = Phi_mu2_eta1*q_acc_eta1;


%% Plot estimation and simulation

figure;
plot(t,adata(11,:),'--r',t,adata_mu2_est(10,:),'k')
set(get(gcf,'CurrentAxes'),'FontName','Times New Roman','FontSize',10)
xlabel('$t$','FontName','Times New Roman','FontSize',10,'Interpreter','latex');
ylabel('$\ddot{d}_{11}(t)$','FontName','Times New Roman','FontSize',10,'Interpreter','latex');
legend('Simulated','Modal expansion','FontName','Times New Roman','FontSize',10,'Interpreter','latex')





%% Functions
function [K,M,freedof]=beam2D(coor,no,fixeddof)
if size(no,1)==1
    dofe=1:6;
else
    dofs=reshape(1:(max(no(:,2))*3),3,[])';
    dofno=num2cell(dofs,2);
    dofe=cell2mat(dofno(no(:,1:2)));
end
ndof=max(max(dofe));
Kglobal = zeros(ndof,ndof);
Mglobal = zeros(ndof,ndof);
for i=1:size(no,1)
    L = sqrt((coor(no(i,2),1)-coor(no(i,1),1))^2 + ...
             (coor(no(i,2),2)-coor(no(i,1),2))^2);
    I = no(i,3);
    E = no(i,4);
    A = no(i,5);
    rho = no(i,6);
    
    Ke = [ E*A/L 0 0 -E*A/L 0 0;
           0 12*E*I/L^3 6*E*I/L^2 0 -12*E*I/L^3 6*E*I/L^2;
           0 6*E*I/L^2 4*E*I/L 0 -6*E*I/L^2 2*E*I/L;
           -E*A/L 0 0 E*A/L 0 0;
           0 -12*E*I/L^3 -6*E*I/L^2 0 12*E*I/L^3 -6*E*I/L^2;
           0 6*E*I/L^2 2*E*I/L 0 -6*E*I/L^2 4*E*I/L ];
    
    Me = rho*A*L/420 * [ 140 0 0 70 0 0;
                         0 156 22*L 0 54 -13*L;
                         0 22*L 4*L^2 0 13*L -3*L^2;
                         70 0 0 140 0 0;
                         0 54 13*L 0 156 -22*L;
                         0 -13*L -3*L^2 0 -22*L 4*L^2 ];
                           
    a = atan((coor(no(i,2),2)-coor(no(i,1),2))/...
            (coor(no(i,2),1)-coor(no(i,1),1)));

    T = [ cos(a) sin(a) 0   0      0     0
          -sin(a) cos(a) 0   0      0     0
          0       0    1   0      0     0 
          0       0    0  cos(a) sin(a) 0
          0       0    0 -sin(a) cos(a) 0
          0       0    0   0      0     1];

    Kel=T'*Ke*T;
    Mel=T'*Me*T;
    Kglobal(dofe(i,:),dofe(i,:))=Kel+Kglobal(dofe(i,:),dofe(i,:));
    Mglobal(dofe(i,:),dofe(i,:))=Mel+Mglobal(dofe(i,:),dofe(i,:));
end
    totdof=1:ndof;
    freedof=totdof;
    freedof(fixeddof)=[];
    K=Kglobal(freedof,freedof);
    M=Mglobal(freedof,freedof);
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



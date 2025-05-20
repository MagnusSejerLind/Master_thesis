clc,clear,close all
set(0,'defaultTextInterpreter','latex');

%% System properties

opt.sysType = "chain";  % ["chain" / "frame"] - Type of system
% opt.method = "TA";      % ["TA"/"ME"] - Virtuel sensing method (Toeplitz's/Modal expansion)
opt.out_type = 0;       % [disp=0 / vel=1 / acc=2] - Define output type
opt.error_mod = 0;      % [0/1] - Include error modeling and noise
opt.nonlinear = 0;      % [0/1] - Include nonlinearties in the system
opt.nonlinType = 1;     % [0=constant / 1=varied] - Define type of nonlineaties
opt.numDOF = 4;          % Number of DOF --ONLY FOR CHAIN SYSTEM
opt

in_dof = [1 3];         % Input DOF
out_dof = [1 3];        % Output DOF
% out_dof = [1 2 3 4 5 6 7 8 9 10 12 15 16 18 19 20 22 23 24];        % Output DOF

%% System modeling

[dof,m,k,xi_int] = systemSetup(opt);
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
% u_mag = 10;
% u = ones(r,N)*u_mag;
% u = u.*sin(t*5);
u = zeros(r,N);
% u(N*0.2) = u_mag;


% Base system
if opt.sysType == "chain"
    if opt.error_mod == 1; [k,m,snr] = modeling_error(k,m); end
    [M,~,K] = chain(m,m*0,k,dof);
end
addBeamError = 1;
if opt.sysType == "frame"; [M,K,dof,snr] = beamStruc(opt,addBeamError); end
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
[alpha,beta] = raylieghDamp(omegaN,xi_int);
C = alpha*M + beta*K;
C_modal = round(Phi'*C*Phi,10);
xi = diag(C_modal) ./ (2*omegaN);

% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,opt.out_type,dt);

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
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);                        % non-linear stiffness force (displacement dependt)

    z_new = Ad*z_old + Bd*u(:,i) - fd_nl - fk_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);
if opt.error_mod == 1
    Y = awgn(Y,snr,'measured');
end


%% FRF est.
% omega=0:2*pi/t(end):2*pi*1/dt;  % frequency vector (time-equliviant)


% U_freq = fft(u(1,:));   % Input frequency domain
% Y_acc_freq=fft(y(1,:));   % displacement at DOF 1 in frequency domain
% h11_est=real(Y_acc_freq)./real(U_freq);  % FRF computation
% semilogy(omega,abs(h11_est))
% xlim([0 100])
% grid

%% FRF analy.

omega = 0:2*pi/t(end):2*pi*1/dt;  % frequency vector (time-equliviant)


% Analytical Receptance frequency response matrix
for j = 1:numel(omega)
    H_an{j} = (-M*omega(j)^2 + C*omega(j)*1i + K)^-1;
end
HH = [H_an{:}];


% Extract h_11
h11_an=HH(1:dof^2:end);

figure()
semilogy(omega,abs(h11_an),'k',LineWidth=2)
xlim([0 50])
currentYLim = ylim;
ylim([currentYLim(1)*0.75  max(real(abs(h11_an)))*2])
title('FRF')
grid minor
xlabel('Hz')



%% chain all dof frf
if opt.sysType == "chain"

    h11 = HH(1:dof^2:end);
    h22 = HH(2:dof^2:end);
    h33 = HH(3:dof^2:end);
    h44 = HH(4:dof^2:end);

    h_diag = [h11; h22; h33; h44];

    figure()
    tl = tiledlayout(2,2)
    title(tl,'Frequency Response Functions, $i=j$','Interpreter','latex')
    for i = 1:dof
        nexttile
        semilogy(omega,abs(h_diag(i,:)),'k',LineWidth=1.5)
        xlim([0 40])
        ylim([1e-5  1e-1])
        grid
        xlabel('$\omega$ [Hz]')
        ylabel(['$|h_{' num2str(i) num2str(i) '}(\omega)|$'])
        yticks([1e-5,1e-4,1e-3,1e-2,1e-1])
    end

end


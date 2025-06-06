clc,clear
close all
set(0,'defaultTextInterpreter','latex');
rng('default')

%% System & options
opt.sysType = "frame";  % ["chain" / "frame"] - Type of system
opt.method = "TA";      % ["TA"/"ME"] - Virtual sensing method (Toeplitz's/Modal expansion)
opt.out_type = 0;       % [0/1/2] - Define output type (0=disp, 1=vel, 2=accel)
opt.numDOF = 8;         % [-int.-] - Number of DOF --ONLY FOR CHAIN SYSTEM
opt.nonlinear = 1;      % [0/1] - Include nonlinearties in the system
opt.nonlinType = 1;     % [0/1] - Define type of nonlineaties (0=constant / 1=varied)
opt.error_mod = 1;      % [0/1] - Include modeling errors and noise
opt.psLoads = 1;        % [0/1] - Apply pseodu loads to model a nonlinear system
opt.condimp = 1;        % [0/1] - Improve condition of Toeplitz's matrix by system truncation
opt.estTempTrunc = 1;   % [0/1] - Temporal truncation of output estimations
opt.plot = 1;           % [0/1] - Plot results
opt.subsetPlot = 0;     % [0/1] - Plots only a subset of the estimated outputs
opt.animate = 1;        % [0/1] - Animates the displacements of the structure
opt.aniSave = 0;        % [0/1] - Save animation
opt.save = 0;           % [0/1] - Save figure
opt

in_dof = [1 2 3];         % Input DOF
out_dof = [3 6 8];        % Output DOF, frame dof: (x,y,θ)

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

% Input 
u_mag = 10;
u = ones(r,N)*u_mag;
u = u.*sin(t*5);
% u = zeros(r,N);
% u(N*0.2:N*0.3) = u_mag;
U = u(:);
if opt.error_mod == 1
    [~,~,snr] = modeling_error(0,0);
    U = awgn(U,snr,'measured');
    u = reshape(U,r,N);
end


% Actucal system
if opt.sysType == "chain"; [M_acc,~,K_acc] = chain(m,m*0,k,dof); end
addBeamError = 0;
if opt.sysType == "frame"; [M_acc,K_acc,dof,snr] = beamStruc(opt,addBeamError); end
[Phi_acc,Lambda_acc] = eig(K_acc,M_acc);    % modal and spectral matrix
[omegaN_acc,i2] = sort(sqrt(diag(Lambda_acc))); % Natural freq. [rad/s]
omegaN_acc = real(omegaN_acc);
Phi_acc = Phi_acc(:,i2);
dd_acc = sqrt(diag(Phi_acc'*M_acc*Phi_acc));
aa_acc = Phi_acc*diag(1./dd_acc);    % Mass-normalized Phi (eigenvec.)
[alpha_acc,beta_acc] = raylieghDamp(omegaN_acc,xi_int);
C_acc = alpha_acc*M_acc + beta_acc*K_acc;
C_modal_acc = round(Phi_acc'*C_acc*Phi_acc,10);

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
Phi = Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi)); % Mass norm M
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
[alpha,beta] = raylieghDamp(omegaN,xi_int);
C = alpha*M + beta*K;
C_modal = round(Phi'*C*Phi,10);
xi = diag(C_modal) ./ (2*omegaN);

% Extended system - full output
in_dof_ex = in_dof;
out_dof_ex = (1:1:dof);
dof_ex = numel(out_dof_ex);
r_ex = numel(in_dof_ex);
ms_ex = numel(out_dof_ex);

% System matricies
[Ad,Bd,Cd,Dd] = systemMatriciesSS_dis(M,K,C,dof,in_dof,out_dof,opt.out_type,dt);
[Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
[Ad_acc,Bd_acc,Cd_acc,Dd_acc] = systemMatriciesSS_dis(M_acc,K_acc,C_acc,dof,in_dof,out_dof_ex,opt.out_type,dt);

% Toeplitz's matricies
[H] = ToeplitzMatrix(N,ms,r,Ad,Bd,Cd,Dd);
[H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

% Nonlinearities
if opt.nonlinear == 1
    nl_mag = 1;
    if opt.nonlinType == 0
        % varied / const. mag.
        cf_nl = nl_mag;    % coeffcient of nonlinear damping
        kf_nl = nl_mag;    % coeffcient of nonlinear stiffness
    else
        cf_nl = rand(1,dof)*nl_mag;
        kf_nl = rand(1,dof)*nl_mag;
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
y = zeros(ms,length(t));
for i = 1:N
    fd_nl(dof+1:end) = cf_nl*z_old(dof+1:end).*abs(z_old(dof+1:end));   % nonlinear damping force (velocity dependent)
    fk_nl(dof+1:end)  = kf_nl*(z_old(1:dof).^3);                        % nonlinear stiffness force (displacement dependent)

    z_new = Ad*z_old + Bd*u(:,i) - fd_nl - fk_nl;
    y(:,i) = Cd*z_old + Dd*u(:,i);
    z_old = z_new;
end
Y = y(:);
if opt.error_mod == 1
% output noise
    Y = awgn(Y,snr,'measured');
    y = reshape(Y,ms,N);
end

% Actual system
z_old_acc = z0;
z_new_acc = zeros(size(z_old_acc));
fd_nl_acc = zeros(size(z_old_acc));
fk_nl_acc = zeros(size(z_old_acc));
y_acc = zeros(dof,length(t));
for i = 1:N
    fd_nl_acc(dof+1:end) = cf_nl*z_old_acc(dof+1:end).*abs(z_old_acc(dof+1:end));   % nonlinear damping force (velocity dependt)
    fk_nl_acc(dof+1:end)  = kf_nl*(z_old_acc(1:dof).^3);                            % nonlinear stiffness force (displacement dependt)

    z_new_acc = Ad_acc*z_old_acc + Bd_acc*u(:,i) - fd_nl_acc - fk_nl_acc;
    y_acc(:,i) = Cd_acc*z_old_acc + Dd_acc*u(:,i);
    z_old_acc = z_new_acc;
end
Y_acc = y_acc(:);

%% Condition improvement

if opt.condimp == 1
% Removes first ms rows, and last r columns in the Toeplitz's matricies

    H = H(ms+1:end, 1:end-r);
    H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);

    U = U(1:end-r);
    u = reshape(U,r,N-1);
    Y = Y(ms+1:end);
end

%% Pseudo loads

if opt.psLoads == 1

    % Full input
    in_dof_FI = 1:1:dof;
    r_FI = numel(in_dof_FI);
    [Ad_FI,Bd_FI,Cd_FI,Dd_FI] = systemMatriciesSS_dis(M,K,C,dof,in_dof_FI,out_dof,opt.out_type,dt);
    [H_FI] = ToeplitzMatrix(N,ms,r_FI,Ad_FI,Bd_FI,Cd_FI,Dd_FI);

    if opt.condimp
        H_FI = H_FI(ms+1:end, 1:end-r_FI);
    end

    Gamma = pinv(H_FI)*(Y - H*U);   % Pseudo loads calc

else
    Gamma = zeros(dof*length(t),1);
    if opt.condimp == 1 && opt.psLoads ~= 1; Gamma = Gamma(1:end-dof); end
end
%% Output estimation

if opt.method == "TA"
% Toeplitz's approach

    % Full input, full output
    in_dof_FIFO = 1:1:dof;
    r_FIFO = numel(in_dof_FIFO);
    ms_FIFO = numel(out_dof_ex);
    [Ad_FIFO,Bd_FIFO,Cd_FIFO,Dd_FIFO] = systemMatriciesSS_dis(M,K,C,dof,in_dof_FIFO,out_dof_ex,opt.out_type,dt);
    [H_FIFO] = ToeplitzMatrix(N,ms_FIFO,r_FIFO,Ad_FIFO,Bd_FIFO,Cd_FIFO,Dd_FIFO);

    if opt.condimp
        H_FIFO = H_FIFO(ms_FIFO+1:end, 1:end-r_FIFO);
    end

    Y_est = H_ex*U + H_FIFO*Gamma;  % output estimations

    if opt.condimp == 1
        y_est = reshape(Y_est, dof, N-1);  
    else
        y_est = reshape(Y_est, dof, N);  
    end
end

if opt.method == "ME"
% Modal expansion

    mu1 = out_dof;              % Observed nodes
    mu2 = 1:dof; mu2(mu1)=[];   % Unobserved nodes
    eta1 = 1:numel(mu1);        % Retained modes

    Phi_mu1_eta1 = aa(mu1,eta1);
    Phi_mu2_eta1 = aa(mu2,eta1);

    Phi_mu1_eta1_PI = (Phi_mu1_eta1'*Phi_mu1_eta1)^-1*Phi_mu1_eta1';    % Pseudo-inverse
    q_out_eta1 = Phi_mu1_eta1_PI*y;

    y_mu2_est = Phi_mu2_eta1*q_out_eta1;    % Estimated output
end

if opt.condimp == 1
    y_acc = y_acc(:,2:end);
    Y_acc = y_acc(:);
    t = t(2:end);
end

%% Temporal truncation
if opt.estTempTrunc == 1
% Truncates the temporal end for the system

    trunc_val = 0.1;   % Percentage of truncation approx
    N_trunc = round(N*(1-trunc_val),-1);    % no. of timesteps for truncated time, rounded to nearest 10

    t = t(1:N_trunc);

    y = y(:,1:N_trunc);
    y_acc = y_acc(:,1:N_trunc);

    if opt.method == "TA"
        y_est = y_est(:,1:N_trunc);
        Y_est = y_est(:);
    elseif opt.method == "ME"
        y_mu2_est = y_mu2_est(:,1:N_trunc);
    end
    Y = y(:);
    Y_acc = y_acc(:);

    u = u(:,1:N_trunc);
    U = u(:);

else
    N_trunc = N

end

%% RMSE

mu1 = out_dof;              % Observed nodes
mu2 = 1:dof; mu2(mu1)=[];   % Unobserved nodes


% Root mean squared error
RMSE = zeros(1,ms);
for i = 1:(dof-ms)
    if opt.method =='TA'; RMSE(i) = (sqrt(mean((y_acc(mu2(i),:) - y_est(mu2(i),:)).^2))); end
    if opt.method == 'ME'; RMSE(i) = sqrt(mean((y_acc(mu2(i),:) - y_mu2_est(i,:)).^2)); end
end
RMSE_tot = mean(RMSE)

if opt.plot == 1
    figure()
    bar(mu2,RMSE,'k')
    xlabel('DOF ')
    ylabel('RMSE')
    title('Root mean squared error')
    grid
    xticks(1:1:dof)
end

% RMSE on dof type
if opt.sysType == "frame"
    x_enum = mod(mu2 - 1, 3) == 0;
    y_enum = mod(mu2 - 1, 3) == 1;
    theta_enum = mod(mu2 - 1, 3) == 2;
    RMSE_dof = [mean(RMSE(x_enum)),mean(RMSE(y_enum)),mean(RMSE(theta_enum))];
    dof_order = ["x","y","$\theta$"];

    if opt.plot == 1
        figure()
        bar(dof_order,RMSE_dof,'k')
        ylabel('RMSE')
        xlabel('DOF type')
        grid
        title('RMSE DOF type')
        ax = gca;
        ax.XAxis.TickLabelInterpreter = 'latex';
    end
end

%% Visualization

if opt.plot == 1
    figure()
    tiledlayout('flow')
    if opt.method == 'TA'; sgtitle("Output estimation - Toeplitz's approach",'Interpreter','latex'); end
    if opt.method == "ME"; sgtitle("Output estimation - Modal expansion",'Interpreter','latex'); end


    if opt.subsetPlot == 1
    % plot only a subset of estimated dof outputs
       
        % subset to plot -- index of mu2 --
        ssmu2 = [1 2 3  21 22 23];

        for i = ssmu2
            nexttile
            plot(t,y_acc(mu2(i),:)','k',LineWidth=2)
            hold on

            if opt.method == 'TA'; plot(t,y_est(mu2(i),:),'r--',LineWidth=2); end
            if opt.method == "ME"; plot(t,y_mu2_est(i,:),'--r',LineWidth=2); end


            legend('Actual output','Estimated output','Interpreter','latex')
            title(sprintf('DOF: %d', mu2(i)));
            grid
            xlabel('Time [s]')
            ylabel(sprintf('Output (%d)', opt.out_type));
            xlim([0 N_trunc*dt])
            ylim([min(y_acc(mu2(i),:))*1.2 max(y_acc(mu2(i),:))*1.5])

        end

    else

        for i = 1:numel(mu2)
            nexttile
            plot(t,y_acc(mu2(i),:)','k',LineWidth=2)
            hold on

            if opt.method == 'TA'; plot(t,y_est(mu2(i),:),'r--',LineWidth=2); end
            if opt.method == "ME"; plot(t,y_mu2_est(i,:),'--r',LineWidth=2); end


            legend('Actual output','Estimated output','Interpreter','latex')
            title(sprintf('DOF: %d', mu2(i)));
            grid
            xlabel('Time [s]')
            ylabel(sprintf('Output (%d)', opt.out_type));
            xlim([0 N_trunc*dt])
            ylim([min(y_acc(mu2(i),:))*1.2 max(y_acc(mu2(i),:))*1.5])
        end
    end
end


% Animate displacements
if opt.animate == 1 && opt.sysType == "frame"
    beamStrucDisplacement(y_est,u,in_dof,opt)
elseif  opt.animate == 1 && opt.sysType == "chain"
    chain_animation(y_est,t,N,opt)
end

%%

if opt.save == 1
    print('Results/NAME_XX', '-dpng');
end

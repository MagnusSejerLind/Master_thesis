clc,clear,
close all
set(0,'defaultTextInterpreter','latex');
rng('default')
opt.plot = 1;           % [0/1]
opt.larm = 1;


opt.condimp = 1;        % [0/1] - Improve Toeplitz's matrix condition by truncation



%% System properties

noOut = 2;  % Number of outputs (sensors)
opt.sysType = "chain";  % ["chain" / "frame"] - Type of system
opt.out_type = 2;       % [disp=0 / vel=1 / acc=2] - Define output type
opt.error_mod = 0;      % [0/1] - Include error modeling and noise
% opt.numDOF = 8;          % Number of DOF --ONLY FOR CHAIN SYSTEM


max_dof = 50;


% For chain dof loop:
optPlace = zeros(max_dof,1);
for  numcodcount= 2:max_dof
% disp(numcodcount/50)
numcodcount
opt.numDOF = numcodcount;


%% System modeling

[dof,m,k,xi] = systemSetup(opt);

% Time
N = 500;
dt = 0.01;
t = 0:dt:(N-1)*dt;

% Base system
if opt.sysType == "chain"
    if opt.error_mod == 1; [k,m,snr] = modeling_error(k,m); end
    [M,~,K] = chain(m,m*0,k,dof);
end
addBeamError = [];
if opt.sysType == "frame"; [M,K,dof,snr] = beamStruc(opt,addBeamError); end
[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
% C_modal = diag(2*xi.*omegaN);
% C = inv(aa)'*C_modal*inv(aa);
[alpha,beta] = raylieghDamp(omegaN,xi);
C = alpha*M + beta*K;
C_modal = round(Phi'*C*Phi,10);

% Extended system ---------not defined same as in main.m
in_dof_ex = (1:1:dof);
r_ex = numel(in_dof_ex);
ms_ex = noOut;


%% chain
if opt.sysType == "chain"


    % m=1
    if noOut == 1
        condNum_calc = NaN(1,max_dof);
        for i = 1:dof

            out_dof_ex = i;

            [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
            [H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);
            if opt.condimp
                H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);
            end
            condNum_calc(i) = cond(H_ex);    % condtion of H_ex:
        end

        [~, minIdx] = min(condNum_calc);
    condNum(dof,:) = condNum_calc(:);

        % fprintf('Optimal sensor placement: DOF %g\n',minIdx)
        optPlace(dof) = minIdx;
    end





    % m=2
    %%% cases med (i,j)=(j,i)' kan skibbes (halver calcu)
    if noOut == 2
        condNum = NaN(dof);
        disp("Completion:")
        for i = 1:dof
            for j = 1:dof
                if i ~= j   % Skip i=j cases

                    out_dof_ex = [i,j];
                    [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
                    [H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);
                    if opt.condimp
                        H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);
                    end
                    condNum(i,j) = cond(H_ex);    % condtion of H_ex:
                end
            end
            disp((i/dof))
        end
        [minValue, linearIdx] = min(condNum(:));
        [rowIdx, colIdx] = ind2sub(size(condNum), linearIdx);
        fprintf('Optimal sensor placement: DOF %g,%g\n',rowIdx,colIdx)

        optPlace(dof,1) = rowIdx;
        optPlace(dof,2) = colIdx;


        if opt.plot == 1
            figure()
            surf(condNum)
            colormap('parula');
            grid minor
            xlim([0 25])
            ylim([0 25])
            ylabel('DOF $i$')
            xlabel('DOF $j$')
            zlabel('$\kappa \left( \tilde{H} \right)$')
            title('Condition number sensor placement - m=2')
            box on
        end
    end



    % m=3
    if noOut == 3
        count =0;
        for i = 1:dof
            for j = 1:dof
                if i ~= j
                    for k = 1:dof
                        if k ~= i
                            if k ~= j
                                count = count+1;
                                out_dof_ex = [i,j,k];
                                [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
                                [H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

                                if opt.condimp
                                    H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);
                                end

                                condNum(numcodcount-2,count) = cond(H_ex);
                            end
                        end
                    end
                end
            end
        end
    end



end


%% Beam frame

if opt.sysType == "frame"

    if noOut == 1

        for i = 1:dof
            out_dof_ex = i
            [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
            [H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

            if opt.condimp
                H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);
            end

            condNum(i) = cond(H_ex);    % condtion of H_ex:
        end
        [~, minIdx] = min(condNum);
        optPlace = minIdx; % OBS. fixed dof not included

        % Determine the node and DOF type
        node = ceil(optPlace / 3);
        dofTypeIndex = mod(optPlace - 1, 3);

        switch dofTypeIndex
            case 0
                dofType = 'x';
            case 1
                dofType = 'y';
            case 2
                dofType = '\theta';
        end
        fprintf('DOF %d is located at node %d and corresponds to the %s degree of freedom.\n', optPlace, node, dofType);

        file_name = 'Results/condNum_m1_beam_disp_N200.mat';
        save(file_name, 'condNum');

        if opt.plot == 1
            figure()
            hold on
            plot(condNum(1:3:end),'.--',LineWidth=1,MarkerSize=30)
            plot(condNum(2:3:end),'.--',LineWidth=1,MarkerSize=30)
            plot(condNum(3:3:end),'.--',LineWidth=1,MarkerSize=30)
            legend('x','y','\theta')
            title('Condition number overs Frame structre DOFs - $m=1$')
            xlabel('Node')
            ylabel('$\kappa \left( \tilde{H} \right)$')
            grid
            ylim([0 6500])
        end
    end



    if noOut == 2
        condNum = NaN(dof);
        disp("Completion:")
        for i = 1:dof
            for j = 1:dof
                if i ~= j   % Skip i=j cases

                    out_dof_ex = [i,j];
                    [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);

                    [H_ex] = ToeplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

                    if opt.condimp
                        H_ex = H_ex(ms_ex+1:end, 1:end-r_ex);
                    end

                    condNum(i,j) = cond(H_ex);    % condtion of H_ex:
                end
            end
            disp((i/dof))
        end
        [minValue, linearIdx] = min(condNum(:));
        [rowIdx, colIdx] = ind2sub(size(condNum), linearIdx);
        fprintf('Optimal sensor placement: DOF %g,%g\n',rowIdx,colIdx)

        optPlace(dof,1) = rowIdx;
        optPlace(dof,2) = colIdx;

        surf(condNum)

        file_name = 'Results/condNum_m2_beam_disp_N200.mat';
        save(file_name, 'condNum');

    end

end

end
%%

if opt.larm == 1
    load gong
    sound(y,Fs)
end


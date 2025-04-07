clc,clear,
close all
set(0,'defaultTextInterpreter','latex');
rng('default')
opt.plot = 1;           % [0/1]
opt.larm = 1;
%% System properties


% optPlace = zeros(50,1);
% for  numcodcount= 1:50
    % disp(numcodcount/50)
    % numcodcount
    % opt.numDOF = numcodcount;

    opt.numDOF = 4;          % Number of DOF
    noOut = 2;  % Number of outputs (sensors)

    opt.sysType = "chain";  % ["chain"] - Type of system
    opt.out_type = 2;       % [disp=0 / vel=1 / acc=2] - Define output type
    opt.error_mod = 0;      % [0/1] - Include error modeling and noise

    %% System modeling

    [dof,m,k,xi] = systemSetup(opt);

    % Time
    N = 500;
    dt = 0.01;
    t = 0:dt:(N-1)*dt;

    % Base system
    [M,~,K] = chain(m,m*0,k,dof);
    [Phi,Lambda] = eig(K,M);    % modal and spectral matrix
    [omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
    omegaN = real(omegaN);
    Phi=Phi(:,i2);
    dd = sqrt(diag(Phi'*M*Phi));
    aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
    C_modal = diag(2*xi.*omegaN);
    C = inv(aa)'*C_modal*inv(aa);


    % Extended system ---------not defined same as in outputEst_general (date: 28/3)
    in_dof_ex = (1:1:dof);
    r_ex = numel(in_dof_ex);
    ms_ex = noOut;


    if noOut == 1
        condNum = NaN(1,dof);
        for i = 1:dof

            out_dof_ex = i;

            [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
            [H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);
            condNum(i) = cond(H_ex);    % condtion of H_ex:
        end

        [~, minIdx] = min(condNum);
        
        % fprintf('Optimal sensor placement: DOF %g\n',minIdx)
        optPlace(dof) = minIdx;
    end





    %%% cases med (i,j)=(j,i)' kan skibbes (halver calcu)

    if noOut == 2
        condNum = NaN(dof);
        disp("Completion:")
        for i = 1:dof
            for j = 1:dof
                if i ~= j   % Skip i=j cases

                    out_dof_ex = [i,j];
                    [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
                    [H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);
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
    end


    


    % end


    %% (og grim)


    % evt. condNum(i+j+k)

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
                                [H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

                                condNum(numcodcount-2,count) = cond(H_ex);
                            end
                        end
                    end
                end
            end
        end
    end
% end


if opt.larm == 1
    load gong
    sound(y,Fs)
end

function [dof,m,k,xi] = systemSetup(opt)

% xi = 0.01;
xi = [0.1, 0.06];   % For first and second mode, to be using in Rayleigh damp. calc. (range: 0:0.1 frame)

if opt.sysType == "chain"
    dof = opt.numDOF;
    m = ones(1,dof)*2;
    k = ones(1,dof)*250;
    % xi = (ones(1,dof)*xi)';
elseif opt.sysType ==  "frame"
    [~, ~, ~, ~, ~, dof, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = loadstructure;
    m = [];
    k = [];
    % xi = (ones(1,dof)*xi)';
else
    disp('System correctly not defined')
end
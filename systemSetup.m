function [dof,m,k,xi] = systemSetup(opt)
% Ensures consistant system constants, mass, stiffness, damping ratio

xi = [0.1, 0.06];   % For first and second mode, to be using in Rayleigh damp. calc. (range: 0:0.1 frame)

if opt.sysType == "chain"
    dof = opt.numDOF;
    m = ones(1,dof)*1;
    k = ones(1,dof)*250;
elseif opt.sysType ==  "frame"
    [~, ~, ~, ~, ~, dof, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = loadstructure;
    m = [];
    k = [];
else
    disp('System correctly not defined')
end
function [dof,m,k,xi] = systemSetup(opt)

xi = 0.01;

if opt.sysType == "chain"
    dof = opt.numDOF;
    m = ones(1,dof)*1;
    k = ones(1,dof)*250;
    xi = (ones(1,dof)*xi)';
elseif opt.sysType ==  "frame"
    [~, ~, ~, ~, ~, dof, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = loadstructure;
    m = [];
    k = [];
    xi = (ones(1,dof)*xi)';
else
    print('System correctly not defined')
end
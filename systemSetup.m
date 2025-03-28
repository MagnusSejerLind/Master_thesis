function [dof,m,k,xi] = systemSetup(opt)

if opt.sysType == "chain"

    dof = opt.numDOF;
    m = ones(1,dof)*1;
    k = ones(1,dof)*250;
    xi = (ones(1,dof)*0.1)';

else 
    print('System not defined')
end
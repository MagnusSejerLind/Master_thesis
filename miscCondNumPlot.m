
clc,clear
close all
set(0,'defaultTextInterpreter','latex')
load('Results\condNum_m1_acc_allConfig_dof2to50.mat')

% condNum(2,:) = NaN;

figure()
surf(condNum)
xlabel('Output DOF location')
ylabel('Number of DOF')
zlabel('$\kappa ( \tilde{H} )$')
grid on
xlim([1 50])
ylim([1 50])
zlim([0 max(max(condNum))])
colormap("jet")
title('Condition number of $\tilde{H}$, $m=1$','Interpreter','latex')
box on

figure()
plot(condNum(50,:),'k',LineWidth=2)
ylabel('$\kappa ( \tilde{H} )$')
xlabel('Output DOF location')
grid
title('Condition number of $\tilde{H}$, $m=1$, $n=50$','Interpreter','latex')
ylim([0 max(condNum(50,:))])







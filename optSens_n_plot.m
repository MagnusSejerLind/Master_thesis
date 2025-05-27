
clc,clear,close all
set(0,'defaultTextInterpreter','latex');

%% m=2
load('Results\optimalSensorPlacement_m2_dof1to25_acc.mat')
dof_vec = 1:1:length(optPlace);

figure()
hold on
plot(dof_vec,dof_vec,'k--',LineWidth=0.5)
grid
ylim([1 25.5])
xlim([1.5 25.5])
xlabel('$n$ DOF')
ylabel('Sensor location DOF')
title('Optimal sensor location, $m=2$')
yticks(1:2:25)
box on
for i = 2:dof_vec(end)

loc=ones(1,2)*i;
line([loc(1), loc(1)], [0, i], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
plot(loc,optPlace(i,:),'xr',LineWidth=1.5,MarkerSize=10)
end

%% m=3

load('Results\optimalSensorPlacement_m3_dof3to12_acc.mat')

dof_vec = 3:1:length(optPlace)+2;


figure()
hold on
plot(dof_vec,dof_vec,'k--',LineWidth=0.5)
grid
ylim([1 12.5])
xlim([2.5 12.5])
xlabel('$n$ DOF')
ylabel('Sensor location DOF')
title('Optimal sensor location, $m=3$')
yticks(1:1:12)
box on
for i =3:12

loc=ones(1,3)*i;
line([loc(1), loc(1)], [0, i], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
plot(loc,cell2mat(optPlace(i-2)),'xr',LineWidth=1.5,MarkerSize=10)
end

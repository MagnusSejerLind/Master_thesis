% ----Phi_acc must be given----
clc,close all

loc_initial = (1:1:dof);
con = zeros(1,dof);

figure('Position', [450, 250, 600, 400]) % [left, bottom, width, height]
tl = tiledlayout('flow');
title(tl,'Chain System Mode Shapes - DOF=4','Interpreter','latex')

for i = 1:dof+1
    nexttile
    hold on

    if i == 1
        plot(loc_initial,con,'.--',MarkerSize=25,LineWidth=1.5,Color=[0.5, 0.5, 0.5])
        text(-0.2, 0, 'Initial state', 'VerticalAlignment', 'middle')
    else
        loc = (1:1:dof)+Phi_acc(:,i-1)';
        plot(loc,con,'k.--',MarkerSize=25,LineWidth=1.5)
    end

    if i ~= 1
        text(-0, 0, ['Mode ', num2str(i-1)], 'VerticalAlignment', 'middle')
    end
    if i == dof+1
        xlabel('Initial DOF location')
    end

    xlim([0.5 4.5])
    grid 
    pbaspect([1 0.1 1]) % X:Y:Z aspect ratio
    yticks([])
    xticks(loc_initial)
end
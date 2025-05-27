% ----Phi_acc must be given----
% also omegaN_acc
clc,close all



freq0 = omegaN_acc / (2 * pi);  % [Hz]

loc_initial = (1:1:dof);
con = zeros(1,dof);

figure('Position', [450, 250, 600, 400]) % [left, bottom, width, height]
tl = tiledlayout('flow');
title(tl,'Mode Shapes - n=4','Interpreter','latex')

for i = 1:dof+1
    nexttile
    hold on

    if i == 1   % initial state
        plot(loc_initial,con,'.--',MarkerSize=25,LineWidth=1.5,Color=[0.5, 0.5, 0.5])
        text(-0.2, 0, 'Initial state', 'VerticalAlignment', 'middle')
        xticks(loc_initial)
    else
        loc = (1:1:dof)+Phi_acc(:,i-1)';
        plot(loc,con,'k.--',MarkerSize=25,LineWidth=1.5)
     text(-0, 0, ['Mode ', num2str(i-1)], 'VerticalAlignment', 'middle')
    xticks(loc_initial)
     xticklabels({})
     
    end

    if i == dof+1
        xlabel('DOF location')
    end

    xlim([0.5 4.5])
    grid 
    pbaspect([1 0.1 1]) % X:Y:Z aspect ratio
    yticks([])
    
end
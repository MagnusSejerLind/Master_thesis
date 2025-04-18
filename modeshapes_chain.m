 % ----Phi_acc must be given----



loc_initial = (1:1:dof)
con = zeros(1,dof);


figure()
tiledlayout('flow')
for i = 1:dof
    nexttile
    loc = (1:1:dof)+Phi_acc(:,i)'
    hold on
    plot(loc_initial,con,'k.--',MarkerSize=20,LineWidth=1)
    plot(loc,con,'.-',MarkerSize=30,LineWidth=2)
    % ylim([-1 2])
    xlim([0.5 4.5])
    grid
    % axis tight
    pbaspect([1 0.1 1]) % X:Y:Z aspect ratio

        title(['Mode ', num2str(i)])

end
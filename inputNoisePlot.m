


figure()
plot(t,u(1:r:end),'k',LineWidth=2)
xlabel('time [s]')
ylabel('u',Rotation=0)
grid
title('Input with noise')
ylim([-1.1 1.1])

figure()
plot(t,u_acc(1:r:end),'k',LineWidth=2)
xlabel('time [s]')
ylabel('u',Rotation=0)
grid
title('Input without noise')
ylim([-1.1 1.1])

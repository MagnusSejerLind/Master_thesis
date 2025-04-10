clc,close all



spacing = 0.05;

Y_norm = Y_acc./max(Y_acc);


hf = figure;
xlim([-spacing spacing*(dof+1)])
ylim([-0.2 0.2])


dof1 = line( ...
    'Parent',gca(), ...
    'XData',[], ...
    'YData',[], ...
    'LineStyle','none', ...
    'Color','r', ...
    'Marker','.', ...
    'MarkerSize',30);

dof2 = line(...
    'Parent', gca(), ...
    'XData', [], ...
    'YData', [], ...
    'LineStyle', 'none', ...
    'Color', 'r', ...
    'Marker', '.', ...
    'MarkerSize', 30);

dof3 = line(...
    'Parent', gca(), ...
    'XData', [], ...
    'YData', [], ...
    'LineStyle', 'none', ...
    'Color', 'r', ...
    'Marker', '.', ...
    'MarkerSize', 30);

dof4 = line(...
    'Parent', gca(), ...
    'XData', [], ...
    'YData', [], ...
    'LineStyle', 'none', ...
    'Color', 'r', ...
    'Marker', '.', ...
    'MarkerSize', 30);


ht = text(0.2,max(Y_acc)*0.5,'','Interpreter','latex');
grid

for i = 1:4:N
    ht.String = ['t = ' num2str(t(i),5) '\ s'];
    set(dof1,'XData',Y_acc(i),'YData',0);
    drawnow 
    set(dof2,'XData',Y_acc(i+1)+spacing,'YData',0);
    drawnow 
    set(dof3,'XData',Y_acc(i+2)+2*spacing,'YData',0);
    drawnow 
    set(dof4,'XData',Y_acc(i+3)+3*spacing,'YData',0);
    drawnow 
end



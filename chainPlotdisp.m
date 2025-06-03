
% Time instance
t_set = 0.2;  % [s] if dt=0.01

%%


if opt.out_type ~= 0; return; end
% disp.


    t_i = find(t == t_set);



    spacing = 1;
    scaling = 100;

    dof = opt.numDOF;


    tiledlayout(2,1)


    for q=1:2
        nexttile

        %%


if q == 1
        y_norm = y_est*scaling;
elseif q == 2
    y_norm = y_acc*scaling;
end

        % figure()
        xlim([-spacing spacing*(dof+1)])
        box on
        axis equal
        ylim([-0.5 0.5])
        yticks([])
        xlim([min(y_norm(1,:))*1.5-spacing/2 (height(y_norm)-1)*spacing+(max(y_norm(end,:))*1.25)])
        xlabel('x')
        title(['Damped mass spring system, scaling factor=', num2str(scaling)]);
        grid

        dofLines = zeros(dof, 1);
        for i = 1:dof
            dofLines(i) = line(...
                'Parent', gca(), ...
                'XData', [], ...
                'YData', [], ...
                'LineStyle', 'none', ...
                'Color', 'r', ...
                'Marker', '.', ...
                'MarkerSize', 30);
        end

        % Initialize line objects for connections between DOFs
        connectionLines = zeros(dof-1, 1);
        for i = 1:dof-1
            connectionLines(i) = line(...
                'Parent', gca(), ...
                'XData', [], ...
                'YData', [], ...
                'LineStyle', '--', ...
                'Color', 'k', ...
                'Marker', 'none', ...
                'LineWidth', 1);
        end

        % Bring DOF markers to the front
        for i = 1:dof
            uistack(dofLines(i), 'top');
        end


        % DOF update
        for j = 1:dof
            set(dofLines(j), 'XData', y_norm(j,t_i)+(j-1)*spacing, 'YData', 0);
        end

        % DOF connestions update
        for j = 1:dof-1
            xData = [y_norm(j,t_i)+(j-1)*spacing, y_norm(j+1,t_i)+j*spacing];
            set(connectionLines(j), 'XData', xData, 'YData', [0 0]);
        end
        drawnow;

    end

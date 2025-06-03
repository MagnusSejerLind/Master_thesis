function chain_animation(y_est,t,N,opt)
% Animates the displacement of the chain system

if opt.out_type == 0

    if opt.aniSave == 1
        rec = VideoWriter('chainDisp_XX.avi');
        open(rec)
    end


    spacing = 1;
    scaling = 50;

    % y_norm = (y_est./max(y_est));
    y_norm = y_est*scaling;

    dof = opt.numDOF;

    figure()
    xlim([-spacing spacing*(dof+1)])
    box on
    axis equal
    ylim([-0.5 0.5])
    yticks([])
    xlim([min(y_norm(1,:))*1.5-spacing/2 (height(y_norm)-1)*spacing+(max(y_norm(end,:))*1.25)])
    xlabel('x')
    title(['Damped mass spring system, scaling factor=', num2str(scaling)]);


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

    % time text
    ht = text(0.2,max(max(y_est))*0.5,'','Interpreter','latex');
    grid

    for i = 1:1:N-1
        ht.String = ['t = ' num2str(t(i), 5) ' s'];
        for j = 1:dof
            set(dofLines(j), 'XData', y_norm(j,i)+(j-1)*spacing, 'YData', 0);
        end

        % Update connection lines
        for j = 1:dof-1
            xData = [y_norm(j,i)+(j-1)*spacing, y_norm(j+1,i)+j*spacing];
            set(connectionLines(j), 'XData', xData, 'YData', [0 0]);
        end

        drawnow;
        if opt.aniSave == 1
            frame = getframe(gcf);
            writeVideo(rec, frame);
        end
        pause(0.01)
    end

    if opt.aniSave == 1
        close(rec)
    end


end
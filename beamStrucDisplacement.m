function beamStrucDisplacement(y,u,opt)
% Plots and animates the displacement of the structure w. input
% y: displacement over time
% u: input over time

if opt.aniSave == 1
    rec = VideoWriter('strucDisp_XX.avi');
    open(rec)
end

ressize = 1;    % timestepping

% scale input magnitude
u_max = max(abs(u(:)));
if u_max > 1
    u = u / u_max;
end


[~, xy, ~, ~, idb, ~, incidenze, l, gamma, ~, ~, ~, ~, posiz, ~, ~] = loadstructure;

figure()
hold on;
scale_factor = 1e2;    % scale for displacements
force_scale = 5e-1;    % scale for force arrows
moment_scale = 2e-1;   % increased scale for rotation arrows
arc_angle = 235;       % degrees of arc sweep for moment

for iStep = 1:ressize:length(y)

    % output and input at time instance
    y_ins = y(:,iStep);
    f_ins = u(:,iStep);

    pause(0.01);
    clf;
    title(['Displacement of structure, scale factor=', num2str(scale_factor)]);
    xlabel('x [m]')
    ylabel('y [m]',Rotation=0)
    xlim([-1 4]);
    ylim([0 3]);
    grid on
    box on

    diseg2_disp_forces(y_ins, f_ins, scale_factor, force_scale, moment_scale, arc_angle, incidenze, l, gamma, posiz, idb, xy);

    if opt.aniSave == 1
        frame = getframe(gcf);
        writeVideo(rec, frame);
    end

end

if opt.aniSave == 1
    close(rec)
end



    function diseg2_disp_forces(y_ins, f_ins, scale_factor, force_scale, moment_scale, arc_angle, incidenze, l, gamma, posiz, idb, xy)
        % Plot deformed shape, external forces and rotation moments per DOF
        n_el = size(incidenze,1);
        n_gdl = length(y_ins);
        n_f = length(f_ins);
        hold on;

        % Plot each beam element
        for k = 1:n_el
            dofIdx = incidenze(k,:);
            xkG = zeros(6,1);
            for j = 1:6
                idx = dofIdx(j);
                if idx>0 && idx <= n_gdl
                    xkG(j) = y_ins(idx);
                end
            end
            xkG = scale_factor * xkG;

            lam = [cos(gamma(k)) sin(gamma(k)) 0; -sin(gamma(k)) cos(gamma(k)) 0; 0 0 1];
            Lambda = [lam zeros(3); zeros(3) lam];
            xkL = Lambda * xkG;

            csi = linspace(0, l(k), 21);
            fu = zeros(6, numel(csi)); fu(1,:) = 1 - csi/l(k); fu(4,:) = csi/l(k);
            u_def = (fu' * xkL)';
            fw = zeros(6, numel(csi));
            fw(2,:) = 2*(csi/l(k)).^3 - 3*(csi/l(k)).^2 + 1;
            fw(3,:) = l(k)*((csi/l(k)).^3 - 2*(csi/l(k)).^2 + csi/l(k));
            fw(5,:) = -2*(csi/l(k)).^3 + 3*(csi/l(k)).^2;
            fw(6,:) = l(k)*((csi/l(k)).^3 - (csi/l(k)).^2);
            w_def = (fw' * xkL)';

            local_coords = [u_def + csi; w_def];
            xyG_elem = lam(1:2,1:2)' * local_coords;
            xy0_elem = lam(1:2,1:2)' * [csi; zeros(1,numel(csi))];

            plot(xy0_elem(1,:)+posiz(k,1), xy0_elem(2,:)+posiz(k,2), 'k--');    % Undeformed i'th element
            plot(xyG_elem(1,:)+posiz(k,1), xyG_elem(2,:)+posiz(k,2), 'k', 'LineWidth',2);   % deformed i'th element
        end

        % Nodal displacements, forces, and moments
        n_nodes = size(idb,1);
        disp_nodes = zeros(n_nodes,2);
        forces = zeros(n_nodes,2);
        moments = zeros(n_nodes,1);
        for k = 1:n_nodes
            idxX = idb(k,1); idxY = idb(k,2); idxR = idb(k,3);
            if idxX>0 && idxX <= n_gdl, disp_nodes(k,1) = y_ins(idxX); end
            if idxY>0 && idxY <= n_gdl, disp_nodes(k,2) = y_ins(idxY); end
            if idxX>0 && idxX <= n_f, forces(k,1) = f_ins(idxX); end
            if idxY>0 && idxY <= n_f, forces(k,2) = f_ins(idxY); end
            if idxR>0 && idxR <= n_f, moments(k)    = f_ins(idxR); end
        end
        disp_nodes = scale_factor * disp_nodes;
        xyG = xy + disp_nodes;

        % plot nodes
        plot(xy(:,1), xy(:,2), 'k.');
        plot(xyG(:,1), xyG(:,2), 'k.', 'LineWidth',2, 'MarkerSize',20);

        % plot force arrows (X and Y)
        fx = force_scale * forces(:,1);
        fy = force_scale * forces(:,2);

        quiver(xyG(:,1), xyG(:,2), fx, zeros(size(fx)), 0, 'r','LineWidth',1.5, 'MaxHeadSize',5);
        quiver(xyG(:,1), xyG(:,2), zeros(size(fy)), fy, 0, 'r', 'LineWidth',1.5, 'MaxHeadSize',5);


        % plot moment as curved arrow per node with larger arrowhead at end
        for k = 1:n_nodes
            M = moments(k);
            if M == 0, continue; end
            center = xyG(k,:);
            r = moment_scale * abs(M);
            ang = linspace(0, sign(M)*arc_angle*pi/180, 100);
            arc_x = center(1) + r*cos(ang);
            arc_y = center(2) + r*sin(ang);
            plot(arc_x, arc_y, 'r', 'LineWidth', 1.5);

            % arrowhead at end of arc
            end_ang = ang(end);
            end_point = [arc_x(end), arc_y(end)];
            tangent = sign(M)*[-sin(end_ang), cos(end_ang)];
            % compute normal to tangent
            normal = [ -tangent(2), tangent(1) ];
            % arrowhead size
            ah_len = r * 0.75;    % Length
            ah_wid = r * 0.6;    % Width
            base = [end_point(1) - ah_len * tangent(1), end_point(2) - ah_len * tangent(2)];
            left = base + ah_wid * normal;
            right = base - ah_wid * normal;
            % draw filled arrowhead triangle
            tip = end_point;                              % pointy tip at arc end
            patch([tip(1), left(1), right(1)], [tip(2), left(2), right(2)], 'r', 'EdgeColor', 'none');
        end

    end
end
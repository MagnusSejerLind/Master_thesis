function beamStrucDisplacement(y,u,in_dof,opt)
% Plots and animates the displacement of the structure w. input
%
% y: displacement over time
% u: input over time
% in_dof: dofs where a input is placed

if opt.out_type == 0

scale_factor = 1e3;    % scale for displacements


    if opt.aniSave == 1
        rec = VideoWriter('strucDisp_XX.avi');
        open(rec)
    end

    ress = 1;    % timestepping resolution

    % scale input magnitude
    u_max = max(abs(u(:)));
    if u_max > 1
        u = u / u_max;
    end

    % all dof input placement
    u_dof = zeros(size(y));
    j=0;
    for i = 1:height(u_dof)
        cont_bi = ismember(i,in_dof);
        if cont_bi == 1
            j = j+1;
            u_dof(i,:) = u(j,:);
        end
    end
    u = u_dof;


    [~, xy, nnod, ~, idb, dof, incidenze, l, gamma, ~, ~, ~, ~, posiz, ~, ~] = loadstructure;

    figure()
    hold on;
    force_scale = 5e-1;    % scale for force arrows
    moment_scale = 2e-1;   % scale for rotation arrows
    arc_angle = 235;       % degrees of arc - momentarrow

    for iStep = 1:ress:length(y)

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

        diseg2_disp_forces(y_ins, f_ins, scale_factor, force_scale, moment_scale, arc_angle, incidenze, l, gamma, posiz, idb, xy, dof, nnod);

        if opt.aniSave == 1
            frame = getframe(gcf);
            writeVideo(rec, frame);
        end

    end

    if opt.aniSave == 1
        close(rec)
    end

end

    function diseg2_disp_forces(y_ins, f_ins, scale_factor, force_scale, moment_scale, arc_angle, incidenze, l, gamma, posiz, idb, xy, dof, nnod)
        % Plot deformed shape, external forces and rotation moments per DOF

        n_el = size(incidenze,1);   % # elements
        n_gdl = dof;                % # dof
        n_f = length(f_ins);        % # inputs
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
            xkG = scale_factor * xkG;   % Global scaled displacement of element k

            lam = [cos(gamma(k)) sin(gamma(k)) 0; -sin(gamma(k)) cos(gamma(k)) 0; 0 0 1];
            Lambda = [lam zeros(3); zeros(3) lam];  % Rotation matrix
            xkL = Lambda * xkG;     % Local scaled displacement of element k

            csi = linspace(0, l(k), 21);

            fu = zeros(6, numel(csi)); fu(1,:) = 1 - csi/l(k); fu(4,:) = csi/l(k);  % Shape function mat. for axial disp.
            u_def = (fu' * xkL)';   % Axial disp. over el. length

            % Shape function mat. for transv. and rot. disp.
            fw = zeros(6, numel(csi));
            fw(2,:) = 2*(csi/l(k)).^3 - 3*(csi/l(k)).^2 + 1;
            fw(3,:) = l(k)*((csi/l(k)).^3 - 2*(csi/l(k)).^2 + csi/l(k));
            fw(5,:) = -2*(csi/l(k)).^3 + 3*(csi/l(k)).^2;
            fw(6,:) = l(k)*((csi/l(k)).^3 - (csi/l(k)).^2);
            w_def = (fw' * xkL)';   % Transv. disp. over el. length

            local_coords = [u_def + csi; w_def];    % local disp.
            xyG_elem = lam(1:2,1:2)' * local_coords;    % deformed k'th element coord.
            xy0_elem = lam(1:2,1:2)' * [csi; zeros(1,numel(csi))];  % undeformed k'th element coord.

            plot(xy0_elem(1,:)+posiz(k,1), xy0_elem(2,:)+posiz(k,2), 'k--');    % Undeformed i'th element
            plot(xyG_elem(1,:)+posiz(k,1), xyG_elem(2,:)+posiz(k,2), 'k', 'LineWidth',2);   % deformed i'th element
        end

        % Nodal displacements, forces, and moments
        n_nodes = nnod;
        disp_nodes = zeros(n_nodes,2);
        forces = zeros(n_nodes,2);
        moments = zeros(n_nodes,1);

        for k = 1:n_nodes
            % Update nodal pos., force, moment

            % x,y,θ dof enumeration
            idxX = idb(k,1);
            idxY = idb(k,2);
            idxR = idb(k,3);

            if idxX <= n_gdl, disp_nodes(k,1) = y_ins(idxX); end    % nod. pos.
            if idxY <= n_gdl, disp_nodes(k,2) = y_ins(idxY); end

            if idxX <= n_f, forces(k,1) = f_ins(idxX); end  % force
            if idxY <= n_f, forces(k,2) = f_ins(idxY); end

            if idxR <= n_f, moments(k) = f_ins(idxR); end   % moment


        end
        disp_nodes = scale_factor * disp_nodes;
        xyG = xy + disp_nodes;

        % plot nodal pos.
        plot(xy(:,1), xy(:,2), 'k.');   % Initial
        plot(xyG(:,1), xyG(:,2), 'k.', 'LineWidth',2, 'MarkerSize',20);

        % plot force arrows (X and Y)
        fx = force_scale * forces(:,1);
        fy = force_scale * forces(:,2);
        quiver(xyG(:,1), xyG(:,2), fx, zeros(size(fx)), 0, 'r','LineWidth',1.5, 'MaxHeadSize',5);
        quiver(xyG(:,1), xyG(:,2), zeros(size(fy)), fy, 0, 'r', 'LineWidth',1.5, 'MaxHeadSize',5);

        % Plot moments
        for k = 1:n_nodes
            if moments(k) == 0, continue; end

            % arc
            center = xyG(k,:);
            r = moment_scale * abs(moments(k));
            ang = linspace(0, sign(moments(k))*arc_angle*pi/180, 100);
            arc_x = center(1) + r*cos(ang);
            arc_y = center(2) + r*sin(ang);
            plot(arc_x, arc_y, 'r', 'LineWidth', 1.5);

            % arrowhead
            ah_len = r * 0.75;    % Length
            ah_wid = r * 0.6;    % Width
            end_ang = ang(end);
            end_point = [arc_x(end), arc_y(end)];
            tangent = sign(moments(k))*[-sin(end_ang), cos(end_ang)];
            normal = [ -tangent(2), tangent(1) ];
            base = [end_point(1) - ah_len * tangent(1), end_point(2) - ah_len * tangent(2)];
            left = base + ah_wid * normal;
            right = base - ah_wid * normal;
            tip = end_point;
            patch([tip(1), left(1), right(1)], [tip(2), left(2), right(2)], 'r', 'EdgeColor', 'none');
        end

    end
end

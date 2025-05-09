clc,close all
clearvars -except y u U
%%

[file_i, xy, nnod, sizee, idb, ndof, incid, l, gamma, m, EA, EJ, T, posiz, nbeam, pr] = loadstructure;



figure()
hold on
xlim([-0.5 3.5])
ylim([0 2.5])

scale_factor = 1E3

for i = 1:1:100

y_ins = y(:,i);
pause(0.01)
clf
xlim([-0.5 3.5])
ylim([0 2.5])


diseg2_disp(y_ins,scale_factor,incid,l,gamma,posiz,idb,xy)


end

function diseg2_disp(y_ins,scale_factor,incidenze,l,gamma,posiz,idb,xy)

% Checking consistency input data
[n_el,~]=size(incidenze);

if length(posiz) ~= n_el
    sprintf('Error: number of nodes in posit matrix differs from number of elements')
    return
end

n_gdl=length(y_ins);

hold on
% Looping on the finite elements
for k=1:n_el
    % Building the nodal displacements vector of each element in the global reference frame
    xkG = zeros(6,1);
    for iri=1:6
        if incidenze(k,iri) <= n_gdl
            xkG(iri,1) = y_ins(incidenze(k,iri));
        end
    end

    % Applying scale factor
    xkG = scale_factor * xkG;

    % Global to Local reference frame rotation
    lambda = [ cos(gamma(k)) sin(gamma(k)) 0. 
              -sin(gamma(k)) cos(gamma(k)) 0.
                     0.      0.     1. ] ;
    Lambda = [ lambda     zeros(3)
              zeros(3)   lambda      ] ;
    xkL = Lambda * xkG;

    % Computing the axial (u) and transversal (w) displacements by means of shape functions
    csi = l(k) * [0:0.05:1];
    fu = zeros(6,length(csi));
    fu(1,:) = 1 - csi / l(k);
    fu(4,:) = csi / l(k);
    u = (fu' * xkL)';

    fw = zeros(6,length(csi));
    fw(2,:) = 2 * (csi / l(k)).^3 - 3 * (csi / l(k)).^2 + 1;
    fw(3,:) = l(k) * ((csi / l(k)).^3 - 2 * (csi / l(k)).^2 + csi / l(k));
    fw(5,:) = -2 * (csi / l(k)).^3 + 3 * (csi / l(k)).^2;
    fw(6,:) = l(k) * ((csi / l(k)).^3 - (csi / l(k)).^2);
    w = (fw' * xkL)';

    % Local to global transformation of the element's deformation
    xyG = lambda(1:2,1:2)' * [u + csi; w];
    undef = lambda(1:2,1:2)' * [csi; zeros(1,length(csi))];

    % Plotting undeformed and deformed element's shape
    plot(undef(1,:) + posiz(k,1), undef(2,:) + posiz(k,2), 'k--')
    plot(xyG(1,:) + posiz(k,1), xyG(2,:) + posiz(k,2), 'k', 'LineWidth', 2)
end

% Looping through the nodes
n_nodi = size(idb,1);
xkG = zeros(n_nodi, 2);
for k=1:n_nodi
    for ixy=1:2
        if idb(k,ixy) <= n_gdl
            xkG(k,ixy) = y_ins(idb(k,ixy));
        end
    end
end
xkG = scale_factor * xkG;
xyG = xkG + xy;

% Plotting original and displaced nodal positions
plot(xy(:,1), xy(:,2), 'k.')
plot(xyG(:,1), xyG(:,2), 'k.', 'LineWidth', 2, 'MarkerSize', 20)

grid on
box on
axis equal
end





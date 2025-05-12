function [file_i, xy, nnod, sizee, idb, dof, incidenze, l, gamma, m, EA, EJ, T, posiz, nbeam, pr] = loadstructure
% imports and interprets the beam structure
% 
% file_i: name of imported file
% xy: location of nodes
% nnod: number of nodes
% sizee: max dimension of structure
% idb: DOF enumeration (global w. fixed dofs included)
% dof: number of DOF
% incidenze: connectivity of DOF and elements
% l: element length [m]
% gamma: element orientation [rad]
% m: element mass [kg]
% EA,EJ,T: material parameters
% posiz: coordinates of first node in elements
% nbeam: number of beams
% pr: property vector of material parameters (if multiple types of beam are used)





% Loads *.inp file
file_i = 'input_struc';
if exist([file_i '.inp'])~=2
    
    return
end
eval(['fid_i = fopen(''', file_i, '.inp'', ''r'');']);

% Reading card NODES
findcard(fid_i, '*NODES');
finewhile = 1;
iconta = 0;

while finewhile == 1
    line = scom(fid_i);
    finewhile = isempty(findstr(line, '*ENDNODES'));

    if finewhile == 1
        tmp = sscanf(line, '%i %i %i %i %f %f')';
        iconta = iconta + 1;

        if iconta ~= tmp(1)
            disp('Errore: nodi non numerati in ordine progressivo');
            break;
        end

        ivinc(iconta, :) = tmp(2:4);
        xy(iconta, :) = tmp(5:6);
    end
end

nnod = iconta;
sizee = sqrt((max(xy(:, 1)) - min(xy(:, 1)))^2 + (max(xy(:, 2)) - min(xy(:, 2)))^2);

% Building IDB matrix
idof = 0;

for i = 1:nnod
    for j = 1:3
        if ivinc(i, j) == 0
            idof = idof + 1;
            idb(i, j) = idof;
        end
    end
end

dof = idof;

for i = 1:nnod
    for j = 1:3
        if ivinc(i, j) == 1
            idof = idof + 1;
            idb(i, j) = idof;
        end
    end
end

% reading card BEAMS
findcard(fid_i, '*PROPERTIES');
finewhile = 1;
iconta = 0;

while finewhile == 1
    line = scom(fid_i);
    finewhile = isempty(findstr(line, '*ENDPROPERTIES'));

    if finewhile == 1
        tmp = sscanf(line, '%i %f %f %f %f')';
        iconta = iconta + 1;
        properties(iconta, :) = tmp(2:5);
    end
end

% reading card BEAMS
findcard(fid_i, '*BEAMS');
finewhile = 1;
iconta = 0;

while finewhile == 1
    line = scom(fid_i);
    finewhile = isempty(findstr(line, '*ENDBEAMS'));

    if finewhile == 1
        tmp = sscanf(line, '%i %i %i %i')';
        iconta = iconta + 1;
        incid(iconta, :) = tmp(2:3);
        pr(1, iconta) = tmp(4);
        m(1, iconta) = properties(tmp(4), 1);
        EA(1, iconta) = properties(tmp(4), 2);
        EJ(1, iconta) = properties(tmp(4), 3);
        T(1, iconta) = properties(tmp(4), 4);
        l(1, iconta) = sqrt((xy(incid(iconta, 2), 1) - xy(incid(iconta, 1), 1))^2 + (xy(incid(iconta, 2), 2) - xy(incid(iconta, 1), 2))^2);
        gamma(1, iconta) = atan2(xy(incid(iconta, 2), 2) - xy(incid(iconta, 1), 2), xy(incid(iconta, 2), 1) - xy(incid(iconta, 1), 1));
        incidenze(iconta, :) = [idb(incid(iconta, 1), :) idb(incid(iconta, 2), :)];
        posiz(iconta, :) = xy(incid(iconta, 1), :);
    end
end

nbeam = iconta;
fclose(fid_i);
end

function [tmp] = findcard(fid, card)
maxiter = 1e5;
frewind(fid);

for i = 1:maxiter
    if feof(fid)
        error(['The following string is not found: ' card]);
    end

    riga = fgets(fid);

    if ~isempty(findstr(riga, card))
        return;
    end
end

error(['The following string is not found: ' card]);
end

function [line] = scom(fid)
line = fgetl(fid);
tmp = sscanf(line, '%s', 1);

while (tmp == '!')
    if feof(fid)
        error('Impossible to find uncommented string');
    end

    line = fgetl(fid);
    tmp = sscanf(line, '%s', 1);
end
end

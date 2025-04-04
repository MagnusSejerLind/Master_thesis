% function

clc,clear,close all
set(0,'defaultTextInterpreter','latex');
opt.plot = 0;
%%


noOut = 1;  % Number of outputs (sensors)
opt.sysType = "beamstruc";  % ["chain"/] - Type of system
opt.out_type = 2;       % [disp=0 / vel=1 / acc=2] - Define output type
opt.error_mod = 0;      % [0/1] - Include error modeling and noise



%% load structure
[file_i, xy, nnod, sizee, idb, dof, incid, l, gamma, m, EA, EJ, T, posiz, nbeam, pr] = loadstructure;

% Draw structure
if opt.plot == 1
    dis_stru(posiz, l, gamma, xy, pr, idb, dof);
    % legend('Principal beam element','','','Reinforcement beam element')
    ylim([-0.25 2.5])
    xlim([-0.5 3.5])
end

% Assemble mass and stiffness matrices
[M, K] = assem(incid, l, m, EA, EJ, T, gamma, idb);

%% Matrix partitioning

MFF = M(1:dof, 1:dof);
MCF = M(dof+1:end, 1:dof);
MFC = M(1:dof, dof+1:end);
MCC = M(dof+1:end, dof+1:end);

KFF = K(1:dof, 1:dof);
KCF = K(dof+1:end, 1:dof);
KFC = K(1:dof, dof+1:end);
KCC = K(dof+1:end, dof+1:end);

M_full = M;
K_full = K;

% Include only free nodes
M = MFF;
K = KFF;

%% System modeling

% Time
N = 500;
dt = 0.01;
t = 0:dt:(N-1)*dt;

xi = (ones(1,dof)*0.01)';


% Obs. assuming classical damping - consider checking validity

[Phi,Lambda] = eig(K,M);    % modal and spectral matrix
[omegaN,i2] = sort(sqrt(diag(Lambda))); % Natural freq.
omegaN = real(omegaN);
Phi=Phi(:,i2);
dd = sqrt(diag(Phi'*M*Phi));
aa = Phi*diag(1./dd);    % Mass-normalized Phi (eigenvec.)
C_modal = diag(2*xi.*omegaN);
C = inv(aa)'*C_modal*inv(aa);


in_dof_ex = (1:1:dof);
r_ex = numel(in_dof_ex);
ms_ex = noOut;


if noOut == 1

    for i = 1:dof
        out_dof_ex = i

        [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
        [H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);

        condNum(i) = cond(H_ex);    % condtion of H_ex:

    end

    [~, minIdx] = min(condNum);
    optPlace = minIdx; % OBS. fixed dof not included


    % Determine the node and DOF type
    node = ceil(optPlace / 3);
    dofTypeIndex = mod(optPlace - 1, 3);

    switch dofTypeIndex
        case 0
            dofType = 'x';
        case 1
            dofType = 'y';
        case 2
            dofType = '\theta';
    end
    fprintf('DOF %d is located at node %d and corresponds to the %s degree of freedom.\n', optPlace, node, dofType);

if opt.plot == 1
    figure()
    hold on
    plot(condNum(1:3:end),'.--',LineWidth=1,MarkerSize=30)
    plot(condNum(2:3:end),'.--',LineWidth=1,MarkerSize=30)
    plot(condNum(3:3:end),'.--',LineWidth=1,MarkerSize=30)
    legend('x','y','\theta')
    title('Condition number of DOFs - Single output')
    xlabel('Node')
    ylabel('Condition number')
    grid 
end



end



if noOut == 2
condNum = NaN(dof);
        disp("Completion:")
        for i = 1:dof
            for j = 1:dof
                if i ~= j   % Skip i=j cases

                    out_dof_ex = [i,j];
                    [Ad_ex,Bd_ex,Cd_ex,Dd_ex] = systemMatriciesSS_dis(M,K,C,dof,in_dof_ex,out_dof_ex,opt.out_type,dt);
                    [H_ex] = TeoplitzMatrix(N,ms_ex,r_ex,Ad_ex,Bd_ex,Cd_ex,Dd_ex);
                    condNum(i,j) = cond(H_ex);    % condtion of H_ex:
                end
            end
            disp((i/dof))
        end
        [minValue, linearIdx] = min(condNum(:));
        [rowIdx, colIdx] = ind2sub(size(condNum), linearIdx);
        fprintf('Optimal sensor placement: DOF %g,%g\n',rowIdx,colIdx)
    
    optPlace(dof,1) = rowIdx;
    optPlace(dof,2) = colIdx;


end








%% Functions

function [file_i,xy,nnod,sizee,idb,dof,incidenze,l,gamma,m,EA,EJ,T,posiz,nbeam,pr]=loadstructure

% Loads *.inp file
% disp(' ');
% file_i=input(' Name of the input file *.inp (without extension) for the structure to be analyzed = ','s') ;
file_i='input_struc';
% disp(' ')
% Check input file
if exist([file_i '.inp'])~=2
    % disp(['File ',file_i,'.inp does not exist' ])
    % file_i=[];
    return
end
nprova=file_i;

% File opening
eval(['fid_i=fopen(''',file_i,'.inp'',''r'');']);

% Reading card NODES
findcard(fid_i,'*NODES')
% line=scom(fid_i)
% finewhile=findstr(line,'*ENDNODES')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDNODES'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i %f %f')';
        iconta=iconta+1;
        if iconta ~=tmp(1)
            disp('Errore: nodi non numerati in ordine progressivo')
            break
        end
        ivinc(iconta,:)=tmp(2:4);
        xy(iconta,:)=tmp(5:6);
    end
end
nnod=iconta;
% disp(['Number structure nodes ',int2str(nnod)])
sizee=sqrt((max(xy(:,1))-min(xy(:,1)))^2+(max(xy(:,2))-min(xy(:,2)))^2);
% End of reading card NODES

% Building IDB matrix
idof=0;
for i=1:nnod
    for j=1:3
        if ivinc(i,j) == 0
            idof=idof+1;
            idb(i,j)=idof;
        end
    end
end
dof=idof;
% disp(['Number of structure d.o.f. ',int2str(dof)])
for i=1:nnod
    for j=1:3
        if ivinc(i,j) == 1
            idof=idof+1;
            idb(i,j)=idof;
        end
    end
end


% reading card BEAMS
findcard(fid_i,'*PROPERTIES')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDPROPERTIES'));
    if finewhile == 1
        tmp=sscanf(line,'%i %f %f %f %f')';
        iconta=iconta+1;
        properties(iconta,:) = tmp(2:5);
    end
end
% disp(['Number of properties ',int2str(iconta)])

% reading card BEAMS
findcard(fid_i,'*BEAMS')
finewhile=1;
iconta=0;
while finewhile ==1
    line=scom(fid_i);
    finewhile=isempty(findstr(line,'*ENDBEAMS'));
    if finewhile == 1
        tmp=sscanf(line,'%i %i %i %i')';
        iconta=iconta+1;
        incid(iconta,:)=tmp(2:3);
        pr(1,iconta) = tmp(4);
        m(1,iconta)  = properties(tmp(4),1);
        EA(1,iconta) = properties(tmp(4),2);
        EJ(1,iconta) = properties(tmp(4),3);
        T(1,iconta) = properties(tmp(4),4);
        l(1,iconta)=sqrt((xy(incid(iconta,2),1)-xy(incid(iconta,1),1))^2+(xy(incid(iconta,2),2)-xy(incid(iconta,1),2))^2);
        gamma(1,iconta)=atan2(xy(incid(iconta,2),2)-xy(incid(iconta,1),2),xy(incid(iconta,2),1)-xy(incid(iconta,1),1));
        incidenze(iconta,:)=[idb(incid(iconta,1),:) idb(incid(iconta,2),:)];
        posiz(iconta,:)=xy(incid(iconta,1),:);
    end
end
nbeam=iconta;
% disp(['Number of beam FE ',int2str(nbeam)])



fclose(fid_i);


end


function tmp=findcard(fid,card)
% Searchimg, within the file specified by fid, and place the file reading
% pointer on the next row. Reteurns 1 if the string is found, 0 otherwise
maxiter=1e5;

frewind(fid);
for i=1:maxiter
    if feof(fid)
        error(['The following string is not found: ' card])
    end
    riga=fgets(fid);
    if ~isempty(findstr(riga,card))
        return
    end
end

error(['The following string is not found: ' card])
end


function line=scom(fid)
line=fgetl(fid);
tmp=sscanf(line,'%s',1);
while(tmp=='!')
    if feof(fid)
        error('Impossible to find uncommented string')
    end
    line=fgetl(fid);
    tmp=sscanf(line,'%s',1);
    %  disp('in scom')
    %  pause
end

end



function dis_stru(posiz,l,gamma,xy,pr,idb,dof)
% Plotting the undeformed structure

xmax = max(xy(:,1));
xmin = min(xy(:,1));
ymax = max(xy(:,2));
ymin = min(xy(:,2));

dx = (xmax - xmin)/100;
dy = (ymax - ymin)/100;
d = sqrt(dx^2 + dy^2);

[PR,~,ic] = unique(pr,'stable');
npr = length(PR);
red     = zeros(1,npr);
green   = cos((0:npr/(npr-1):npr)*pi/2/npr).^2;
blue    = cos(((0:npr/(npr-1):npr)-npr)*pi/2/npr).^2;
colori = [red' green' blue'];
colori = colori(ic,:);

figure();
hold on;
% Step 2: elements
for i=1:length(posiz)
    xin=posiz(i,1);
    yin=posiz(i,2);
    xfi=posiz(i,1)+l(i)*cos(gamma(i));
    yfi=posiz(i,2)+l(i)*sin(gamma(i));
    colore = colori(i,:);
    plot([xin xfi],[yin yfi],'linewidth',2,'color','k');
    %plot([xin xfi],[yin yfi],'b','linewidth',2);
end
grid on; box on;

% Step 1: nodal positions
plot(xy(:,1),xy(:,2),'r.','markersize',20);

triangolo_h = [ 0 0; -sqrt(3)/2 .5; -sqrt(3)/2 -.5; 0 0]*d*2;
triangolo_v = [ 0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2;
triangolo_r = [0 0; .5 -sqrt(3)/2; -.5 -sqrt(3)/2; 0 0]*d*2 * [sqrt(2)/2 -sqrt(2/2); -sqrt(2)/2 -sqrt(2)/2];


hold on
for ii = 1:size(xy,1)
    %rectangle('Position',[xy(ii,1)-d/2 xy(ii,2)-d/2 d d],'curvature',1,'edgecolor','r','linewidth',3);
    text(xy(ii,1) + d, xy(ii,2) + d,num2str(ii));
    if (idb(ii,1) > dof)
        fill(xy(ii,1) + triangolo_h(:,1),xy(ii,2) + triangolo_h(:,2),'b');
    end
    if (idb(ii,2) > dof)
        fill(xy(ii,1) + triangolo_v(:,1),xy(ii,2) + triangolo_v(:,2),'b');
    end
    if (idb(ii,3) > dof)
        fill(xy(ii,1) + triangolo_r(:,1),xy(ii,2) + triangolo_r(:,2),'b');
    end
end

axis equal


title('Beam Structure')
end


function [M,K] = assem(incidenze,l,m,EA,EJ,T,gamma,idb)

% Checking consistency input data
[n_el,nc]=size(incidenze);
if nc ~= 6
    disp('Error: number of columns of incidence matrix different from 6')
    return
end
if length(l) ~= n_el
    sprintf('Error: number of elements in l different from n')
    return
end
if length(m) ~= n_el
    sprintf('Error: number of elements in m different from number of elements')
    return
end
if length(EA) ~= n_el
    sprintf('Error: number of elements in EA different number of elements')
    return
end
if length(EJ) ~= n_el
    sprintf('Error: number of elements in EJ differenc number of elements')
    return
end
if length(T) ~= n_el
    sprintf('Error: number of elements in T different number of elements')
    return
end
if length(gamma) ~= n_el
    sprintf('Error: number of elements in alpha different number of elements')
    return
end

if min(min(incidenze)) ~= 1
    sprintf('Error: dof sequence does not start from 1')
    return
end

% Limiting the total dof of the structure
n_dof=max(max(idb));
% if n_dof > 100
%     sprintf('Errore: dof > 100, not allowed!')
%     return
% end

% Assembling matrices M and K
M=zeros(n_dof,n_dof);
K=zeros(n_dof,n_dof);
for k=1:n_el
    [mG,kG] = el_tra(l(k),m(k),EA(k),EJ(k),T(k),gamma(k));
    for iri=1:6
        for ico=1:6
            i1=incidenze(k,iri);
            i2=incidenze(k,ico);
            M(i1,i2)=M(i1,i2)+mG(iri,ico);
            K(i1,i2)=K(i1,i2)+kG(iri,ico);
        end
    end
end

end

function [mG,kG] = el_tra (l,m,EA,EJ,T,gamma)
% mass matrix in the local reference frame
mL = m*l*[1./3.  0.          0.          1./6.  0.          0.
    0.     13./35.     11.*l/210.  0.     9./70.      -13*l/420.
    0.     11.*l/210.  l^2/105.    0.     13*l/420.   -l^2/140.
    1./6.  0.          0.          1./3.  0.          0.
    0.     9./70.      13*l/420.   0.     13./35.     -11.*l/210.
    0.     -13*l/420.  -l^2/140.   0.     -11.*l/210. l^2/105.   ] ;

% stiffness matrix in the local reference frame
% contribution due to the axial deformation
kL_ax = EA/l* [ 1 0 0 -1 0 0
    0 0 0  0 0 0
    0 0 0  0 0 0
    -1 0 0  1 0 0
    0 0 0  0 0 0
    0 0 0  0 0 0 ] ;

% contribution due to bending deformation
kL_fl = EJ * [ 0.    0.       0.      0.    0.       0.
    0.  +12./l^3+T./l./EJ  6./l^2   0.  -12./l^3-T./l./EJ  6./l^2
    0.   6./l^2  +4./l     0.   -6./l^2  +2./l
    0.    0.       0.      0.    0.       0.
    0.  -12./l^3-T./l./EJ  -6./l^2  0.  +12./l^3+T./l./EJ  -6./l^2
    0.   6./l^2  +2./l    0.   -6./l^2  +4./l    ] ;

kL = kL_ax+kL_fl ;

% matrix transformation from local to global reference frame
% rotation matrix 3x3 (single node DoF)
lambda = [ cos(gamma) sin(gamma) 0.
    -sin(gamma) cos(gamma) 0.
    0.         0.      1.] ;

% rotation matrix 6x6 (element DoF)
Lambda = [ lambda     zeros(3)
    zeros(3)   lambda      ] ;

mG = Lambda' * mL * Lambda ;
kG = Lambda' * kL * Lambda ;

end

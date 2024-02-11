function [Kff,Mff,EG,PhiSE,Wnff,NS]=getEig_RitzSE(ipoltype,pDeg,BC_IDS,Ku,Mu,Ainv,xMesh,yMesh,plot_flag)
%% Solve eigenvalue problem for nodal representation of Ritz super-element

% Evaluation mesh
xMesh2=reshape(xMesh,1,numel(xMesh));
yMesh2=reshape(xMesh,1,numel(yMesh));
[W,dWx,dW,dWxy,dWxx,dWyy]=FdF(ipoltype,pDeg,xMesh2,yMesh2);

% Imposing Boundary conditions - Ritz Super Element
nDOF=(pDeg+1)*(pDeg+1);
IDF=1:nDOF      % IDs for all DOF
IDF(BC_IDS')=[]; % Empty matrix for constraied DOF

Kff = Ku(IDF,IDF);%  --------- (32A)')
Mff = Mu(IDF,IDF);%  --------- (32B)')

Kff = double(Kff);
Mff = double(Mff);

KM_matrix_qualities(Kff,Mff,plot_flag)

[QQ,EE]=eig(Kff,Mff);
[EG,ii]=sort(diag(EE));
Q=QQ(:,ii);
Wnff=sqrt(EG)

%% Recover shape functions - Ritz SE
disp(['Size of [Q]=', num2str(size(Q))])
disp(['Size of [W]=', num2str(size(W))])
disp(['Size of [Ainv]=', num2str(size(Ainv))])
disp(['Size of [xMesh]=',num2str(size(xMesh))])
disp('Recover Shape function using NS=W*Ainv')

%NS=W'*Ainv;
NS=W*Ainv;
NS=NS(:,IDF);
disp(['Size of [NS]=',num2str(size(NS))])

%% Recover eigenfunctions - Ritz SE
for im=1:3
    qm=Q(:,im);
    PhiSE(:,im)=NS*qm;
    disp(['Size of [PhiSE]=',num2str(size(PhiSE))])
    if(plot_flag>=1)
        % Plot 2-D modes shapes
        %WW=reshape(PhiSE,size(xMeshSE));
        WW=reshape(PhiSE(:,im),size(xMesh));
        %disp('size of WW='); size(WW)
        figure
        %surfc(xMeshSE,yMeshSE,WW)
        surfc(xMesh,yMesh,WW)
        shading interp
        colorbar
    end
end
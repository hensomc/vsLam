function [EgN,Ks,Ms,PHI,EG]=beamSolve(ipol,pDeg,intMethod,bcMethod,plot_flag,L,b,h,E,rho,k_spring)
%beamSolve - main driver to solve beam vibration & buckling
%
% Inputs:
% pDeg      = degree for displacement shape function
% plot_flag = show or hide graphics (boolean true/false)
% 
% Outputs:
% EgN = vector of eigenvalues
% Ks  = reduced stiffness matrix
% Ms  = reduced mass matrix

%% Beam geometry & properties
EI=E*b*(h^3)/12.;
rhoA=rho*b*h/386.4;

%% Generate mesh of grids (For defining BCs)
% nGridX = pDeg+1;
% xDomain=linspace(0,L,nGridX);
% xMesh=xDomain;
% gridID=[ (1:nGridX)' xMesh(:) ];
nGridS = pDeg+1;
sDomain=linspace(-1,1,nGridS);
sMesh=sDomain;
sGridID=[ (1:nGridS)' sMesh(:) ]

%% Grid IDs for BCs applied to Beam EndPts

% Edge Grid IDs 
edge1GID=[1];
edge3GID=[pDeg+1];

% One end constrained
BC=[edge1GID]

% Both ends constrained
%BC=[edge1GID edge3GID]

% Define XY coords of BC points
sBCXY=[sGridID(BC,2)]


%% Define Integration points

% Uniform mesh - Simpson rule
nISubDiv = 50;                                    % number of integration subdivisions
xIDomain = linspace(0,L,nISubDiv);                % integration domain x-dir
sIDomain = xIDomain*(2/L) - 1
%dx = xIDomain(2)-xIDomain(1)                      % dx = distance between 2 integration points
ds = (sIDomain(2)-sIDomain(1))                    % ds = distance between 2 integration points
xIMesh = xIDomain;
sIMesh = sIDomain;

% Define  integration points
% IID = 1:nISubDiv-1;            % integration point IDs
% xI  = xIDomain(:) + dx/2;      % x-midpoint of integration element
%IID  = find(xI<L);             % filter out points on boundary x = L
%xI   = xI(IID);
sI  = sIDomain(:) + ds/2;
sIID = find(sI<1);             % filter out points on boundary x = L
sI   = sI(sIID);



%% Integration pts & wts - Gauss quadrature
if( intMethod == 2 )
    Ng=pDeg;  % #gauss pts

    % Gauss pts and weights
    [gp,gw]=lgwt(Ng,-1,1);
    sI=gp;
end

%% Stiffness and Mass matrices
Kc=0;
Mc=0;
wt=1.0;
for i=1:length(sI)
    if(intMethod==2)
        wt=gw(i);
    end
    % Evaluate F and its derivatives at integration point
    [F,dFx,dFxx]=FdF_1D(ipol,pDeg,sI(i));
    
    % Stiffness matrix
    Kc=Kc + dFxx'*EI*dFxx*(4/L^2)^2*(L/2)*wt;
    
    % Mass matrix
    Mc=Mc + F'*rhoA*F*(L/2)*wt;
end
Kc
Mc
%% Check symmetry of Stiffness & Mass matrices - based on term with max diff
maxKERR=max(max(abs(Kc-Kc')))
maxMERR=max(max(abs(Mc-Mc')))

disp('Size of K = '); size(Kc)
disp('Rank of K = '); rank(Kc)
disp('Size of M = '); size(Mc)
disp('Rank of M = '); rank(Mc)
format long
disp('K matrix eigenvalues'); eig(Kc)
disp('M matrix eigenvalues'); eig(Mc)

%% BCs using Null Space method
if(bcMethod == 1)
    IBC=0;
    for ibc=1:length(BC)
        sn = sBCXY(ibc,1);     
        [F,dFx,dFxx]=FdF_1D(ipol,pDeg,sn);
        
        IBC=IBC+1; %  Note: IBC and ibc are 2 variables
        F1(IBC,:)=F;  % constrain w=0
        
        IBC=IBC+1;
        F1(IBC,:)=dFx;  % constrain dw/dx=0
        
    end
    FBC=F1;
    T=null(FBC);
    Ks=double(T'*Kc*T);
    Ms=double(T'*Mc*T)
end

%% Spring stiffnesses for edge BCs
if(bcMethod == 2)
    Kbc=0;
    for ibc=1:length(BC)
        sn = sBCXY(ibc,1);
        [F,dFx,dFxx]=FdF_1D(ipol,pDeg,sn);

        %disp('BCs: size of F='); size(F)
        
        % transverse z-dir springs
        FBC1 = double(F);  % using FdF alt method
        Kbc = Kbc + FBC1'*FBC1;
        
        % torsional springs dw/dx
        FBC2=double(dFx);  %using FdF alt method
        Kbc = Kbc + FBC2'*FBC2;
        
    end
    Kbc = k_spring*Kbc
    disp('Size of Kspr = '); size(Kbc)
    Ks=Kc + Kbc;
    Ms=Mc;
end

%% Check condition numbers
Kcond=cond(Ks)
Mcond=cond(Ms)
size(Ks)
size(Ms)

%% Solve eigenvalue problem
[PHI,EG]=eig(Ks,Ms);
EG = EG;   

%% Sort eigensolutions    
[Eg,Is]=sort(diag(EG));
Phi=PHI(:,Is);
    

%% Plot Modes
NMode=min([7 length(Eg)]);
cc=hsv(NMode);
for im=1:NMode
    
    % Recover polynomial coefficienta
    if(bcMethod==1)
        C=T*Phi(:,im);
    elseif(bcMethod==2)
        C=Phi(:,im);
    end
    
    % Evalaute mode shape at output points (use integration pt coords)
    [FOUT]=FdF_1D(ipol, pDeg, sIMesh(:));
    [FOUTx]=FdF_1D(ipol, pDeg, xIMesh(:));
    
    Wmode=FOUT*C;
    if(im==1)
        figure,plot(xIMesh(:),Wmode(:),'color',cc(im,:));hold;
    else
        plot(xIMesh(:),Wmode(:),'color',cc(im,:));
    end
end
title(strcat('Ritz Beam Free Vibration Mode Shapes: ', 'P = ', int2str(pDeg)));
xlabel('x-position(in)'),ylabel('Transverse Deflection');
legend('Mode1', 'Mode2', 'Mode3', 'Mode4', 'Mode5', 'Mode6', 'Mode7');


%% Compute errors from ref values
%eigErrBeam(EgN, L, EI, rhoA);
%% Recover shape funtions
syms x
N=pDeg+1; II=0:(N-1);
Psi=x.^II
dPsi=diff(Psi,x)
dPsi2=diff(dPsi,x)

D=EI;y=dPsi2;
Kc_e=double(int(y.'*D*y,x,0,L)) %   -----(6)')
Mc_e=double(int(Psi.'*rhoA*Psi,x,0,L)) %  ------(3)')

% diffK=Kc-Kc_e
% diffM=Mc-Mc_e

%% Ritz-Super-Elament Method
XN=linspace(0,1,(pDeg+1)/2)*L;
%XN=linspace(0,1,(pDeg+1)/2);
clear A
id=0;
for i=1:(pDeg+1)/2
    Xi=XN(i);
    %[F,dFx,dFxx]=FdF_1D(ipol,pDeg,Xi);
    id=id+1;
    A(id,:)=subs(Psi,x,Xi);
    %A(id,:)=F;
    id=id+1;
    %A(id,:)=dFx;
    A(id,:)=subs(dPsi,x,Xi);
end
A
t=inv(A)
Mc
Mu=t'*Mc*t
Ku=t'*Kc*t


%% (d) Impose Nodal Boundary conditions
IDS=1:2;
IDF=1:pDeg+1;
IDF(IDS)=[]

Kff = Ku(IDF,IDF);%  --------- (32A)')
Mff = Mu(IDF,IDF) ;% --------- (32B)')

[QQ,EE]=eig(double(Kff),double(Mff))
[EG,ii]=sort(diag(EE));
Q=QQ(:,ii);
% Echo eigenvalues
WnS=sqrt(EG)  % Super Ritz
EgN=sqrt(Eg)  % Ritz


%xn=linspace(0,L,101);
xn=xIMesh;
Psin(1,:)=xn.^0;
for i=2:N
    Psin(i,:)=double(subs(Psi(i),x,xn));
end
% for i=1:pDeg+1
%     [F,dFx,dFxx]=FdF_1D(ipol,pDeg,xIMesh);
%     Psin(i,:)=F;
% end
NS=Psin'*t;   % (A)
NS=NS(:,IDF); % (B)

%% Recover eigenfunctions
for im=1:NMode
    qm=Q(:,im);
    Phix=NS*qm;
    if(im==1)
        figure,plot(xIMesh(:),Phix(:),'color',cc(im,:));hold;
    else
        plot(xIMesh(:),Phix(:),'color',cc(im,:));
    end
end
title(strcat('Super Ritz Beam Free Vibration Mode Shapes: ', 'P = ', int2str(pDeg)));
xlabel('x-position(in)'),ylabel('Transverse Deflection');
legend('Mode1', 'Mode2', 'Mode3', 'Mode4', 'Mode5', 'Mode6', 'Mode7');

end



%%
function eigErrBeam = eigErrBeam(EgN,L,EI,rhoA)
%chk_eigErr(EgN)
%    EgN = Input eigenvalues
%    eigB = Baseline eigenvalues
%    ilam = analysis case to be checked (see main driver)

% Copy input eigenvalues to match baseline array size
length( EgN )
if( length(EgN) >= 5 )
  eigCopy(:,1) = EgN(1:5)
else
  eigCopy(:,1) = EgN
end

% Get baseline eignvalues
lambda=[1.875 4.694 7.855 10.995 14.137];
for (i=1:length(eigCopy() ))
  eigB(i,1) = ((lambda(i)^2)*sqrt(EI/rhoA))/(L^2);
end
eigB

eigErr = abs((eigCopy - eigB)./eigB)*100;
figure,scatter( eigB, eigCopy);

title( {'\bfCalculated vs Baseline Eigenvalues'; get_prog_vers } )
xlabel('\bfBaseline');
ylabel('\bfCalculated');
grid on;

% set x&y axes equal
max_eig = max([max(eigB) max(real(eigCopy))]);
min_eig = min([min(eigB) min(real(eigCopy))]);
axis([0, max_eig, 0, max_eig]);

% reference line
refline(1,0);

end
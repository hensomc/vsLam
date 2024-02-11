%function PlateVibSymCheck(panXY,ilam,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType)
%function [EgN_ret,Ks_ret,Ms_ret,PHI_ret,EG_ret] = PlateVibSymCheck(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType)
function [EgN_ret,Ks_ret,Ms_ret,PHI_ret,EG_ret] = PlateVibSymCheck(panXY,lamDATA,solDATA)

% Symbolic with exact integration
%[EgN_e,Ks_e,Ms_e,PHI_e,EG_e]=PlateVibSym(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType);
[EgN_e,Ks_e,Ms_e,PHI_e,EG_e]=PlateVibSym(panXY,lamDATA,solDATA);

% Symbolic with numerical integration and scaled     
%[EgN_i1,Ks_i1,Ms_i1,PHI_i1,EG_i1]=PlateVibSymIso1(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType);
%[EgN_i1,Ks_i1,Ms_i1,PHI_i1,EG_i1]=PlateVibSymIso1(panXY,lamDATA,ipoltype,solDATA);

% Symbolic with numerical integration and {x,y}       
%[EgN_i2,Ks_i2,Ms_i2,PHI_i2,EG_i2]=PlateVibSymIso2(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType);
%[EgN_i2,Ks_i2,Ms_i2,PHI_i2,EG_i2]=PlateVibSymIso2(panXY,lamDATA,solDATA);

%Wn_compare=[sqrt(EgN_e) sqrt(EgN_i1) sqrt(EgN_i2)]

EgN_ret=EgN_e; Ks_ret=Ks_e; Ms_ret=Ms_e; PHI_ret = PHI_e; EG_ret = EG_e;
%EgN_ret=EgN_i1; Ks_ret=Ks_i1; Ms_ret=Ms_i1; PHI_ret = PHI_i1; EG_ret = EG_i1;
%EgN_ret=EgN_i2; Ks_ret=Ks_i2; Ms_ret=Ms_i2; PHI_ret = PHI_i2; EG_ret = EG_i2;



%function [ EgN,Ks,Ms,PHI,EG ] = PlateVibSym(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType)
function [ EgN,Ks,Ms,PHI,EG ] = PlateVibSym(panXY,lamDATA,solDATA)
%PlateVibSym: Solves plate vibration using symbolic notation
%   Detailed explanation goes here

display('============== Solve in physical coords {x,y} ===============');

[ipoltype,pDeg,Mg,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

% Symbolic variables
syms x y
syms a b th       % plate length, width, thickness
syms rho nu D   % plate properties (D=Et^3/(1-nu^2))

% Symbolic polynomial basis function
[F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_Sym(ipoltype,pDeg); %-----(1)

% Compute [Mc] and [Kc] in (x,y) coords
Mc=int(F.'*rho*th*F,x,[0,a]);  %-----(2)
Mc=int(Mc,y,[0,b])

FxxFyy_sym = dFxx'*dFyy + dFyy'*dFxx;
%Kc=int(D*(dFxx.'*dFxx + dFyy.'*dFyy + 2*nu*dFxx.'*dFyy + 2*(1-nu)*dFxy.'*dFxy),x,[0,a]);  %-----(3)
Kc=int(D*(dFxx.'*dFxx + dFyy.'*dFyy + 2*nu*FxxFyy_sym + 2*(1-nu)*dFxy.'*dFxy),x,[0,a]);  %-----(3)
Kc=int(Kc,y,[0,b])

% RitzSE
if(RitzSE==1)
    id=0;
    
    if(ipoltype == 2)
        nNodeSE=(pDeg+1)*(pDeg+1)-12 + 4;
    elseif(ipoltype == 1 || ipoltype == 3)
        nNodeSE=(pDeg)*(pDeg)-12 + 4 - 1;
    end
    for i=1:nNodeSE
        id=id+1;
        A(id,:)=F;  % All nodes
        if(i<=4)    % Only corner nodes get slope derivatives
            id=id+1;
            A(id,:)=dFx;
            id=id+1;
            A(id,:)=dFy;
        end
    end
    nNodeSE
    A

    psiSE=[-1  1  1 -1  0  1  0 -1];
    etaSE=[-1 -1  1  1 -1  0  1  0];
    
    id=0;
    for i=1:nNodeSE
        id=id+1;
        A(id,:)=subs(A(id,:),[x y],[psiSE(i) etaSE(i)]);
        if(i<=4)    % Only corner nodes get slope derivatives
            id=id+1;
            A(id,:)=subs(A(id,:),[x y],[psiSE(i) etaSE(i)]);
            id=id+1;
            A(id,:)=subs(A(id,:),[x y],[psiSE(i) etaSE(i)]);
        end
    end
    A
    rank(A)
    eig(A)
    Ainv=inv(A)
end  % End Ritz SE


% Numerical [Mc] and [Kc]
[a_n, b_n, t_n, rho_n, D_n, nu_n, abd]=get_numeric_vals(lamDATA,panXY);
Mc=double(subs(Mc, [rho,th,a,b], [rho_n,t_n,a_n,b_n]))  %-----(4)
Kc=double(subs(Kc, [D,nu,a,b], [D_n,nu_n,a_n,b_n]))  %-----(5)

% BCs - clamped at x=0
numBC=pDeg+1;
yBC=linspace(0,b_n,numBC);
xBC=linspace(0,0,numBC);
IBC=0
for i=1:numBC
    IBC=IBC+1;
    FBC(IBC,:)=double(subs(F,[x y],[xBC(i) yBC(i)]));
    IBC=IBC+1;
    FBC(IBC,:)=double(subs(dFx,[x y],[xBC(i) yBC(i)]));
end
T=null(FBC);  %-----(6)
Ks=T'*Kc*T  %-----(7)
Ms=T'*Mc*T  %-----(8)

% Eigenvalue solution
[PHI,EG]=eig(Ks,Ms);  %-----(9)
EgN=sort(diag(EG));
Wn=sqrt(EgN)
%eigError=eigErr(Wn,lamDATA.ilam,bcType,a_n,b_n,abd,ipoltype,pDeg)
eigError=eigErr(panXY,solDATA,lamDATA,Wn);



%%
%function [ EgN,Ks,Ms,PHI,EG ] = PlateVibSymIso1(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType)
function [ EgN,Ks,Ms,PHI,EG ] = PlateVibSymIso1(panXY,lamDATA,solDATA)
%PlateVibSymIso: Solves plate vibration using symbolic notation -
%                in isoparametric coordinates
%   Detailed explanation goes here

display('============== Solve in normalized coords {s,t} ===============');

[ipoltype,pDeg,Mg,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

% Symbolic variables
syms x y
syms s t
syms a b th       % plate length, width, thickness
syms rho nu D   % plate properties (D=Et^3/(1-nu^2))

% Symbolic polynomial basis function
[F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_Sym(ipoltype,pDeg) %-----(1)
F=subs(F, [x y], [s t]);
dFs=subs(dFx, [x y], [s t]);
dFt=subs(dFy, [x y], [s t]);
dFss=subs(dFxx, [x y], [s t]);
dFtt=subs(dFyy, [x y], [s t]);
dFst=subs(dFxy, [x y], [s t]);

% Gauss Legendre points and weights
if(intMethod == 2)
    order_1d=[pDeg+1,pDeg+1];
    order_nd = prod ( order_1d(1:2) );
    [gl_pts,gl_wts]=gl_grid( 2, order_1d, order_nd );
    if(plot_flag ==3)
        gl_grid_display2( order_1d, gl_pts );
    end
    sI=gl_pts(1,:);
    tI=gl_pts(2,:);
end

% Compute [Mc] and [Kc] in (s,t) coords
Kc=0;Mc=0;
for i=1:length(sI);
    psi=sI(i); eta=tI(i); gl_wt=gl_wts(i);
    
    H=subs(F,[s,t],[psi,eta]); Fs=subs(dFs,[s,t],[psi,eta]); Ft=subs(dFt,[s,t],[psi,eta]);
    Fst=subs(dFst,[s,t],[psi,eta]); Fss=subs(dFss,[s,t],[psi,eta]); Ftt=subs(dFtt,[s,t],[psi,eta]);
    FssFtt_sym = Fss'*Ftt + Ftt'*Fss;

    J = jacob2D_Iso( psi, eta, panXY );

%     Kc = Kc + D*(Fss.'*Fss + Ftt.'*Ftt + 2*nu*FssFtt_sym + 2*(1-nu)*Fst.'*Fst)*det(J)*gl_wt;
%     Mc = Mc + H.'*rho*th*H*det(J)*gl_wt;
    Kc = Kc + D*(Fss.'*Fss*16/a^4 + Ftt.'*Ftt*16/b^4 + (2*nu*FssFtt_sym + 2*(1-nu)*Fst.'*Fst)*16/(a^2*b^2))*det(J)*gl_wt;
    Mc = Mc + H.'*rho*th*H*det(J)*gl_wt;
end

% Get numerical values
[a_n, b_n, t_n, rho_n, D_n, nu_n, abd]=get_numeric_vals(lamDATA,panXY);
Mc=double(subs(Mc, [rho,th,a,b], [rho_n,t_n,a_n,b_n]))  %-----(4)
Kc=double(subs(Kc, [D,nu,a,b], [D_n,nu_n,a_n,b_n]))  %-----(5)

% BCs - clamped at x=0
numBC=pDeg+1;
yBC=linspace(-1,1,numBC);
xBC=linspace(-1,-1,numBC);
IBC=0
for i=1:numBC
    IBC=IBC+1;
    FBC(IBC,:)=double(subs(F,[s t],[xBC(i) yBC(i)]));
    IBC=IBC+1;
    FBC(IBC,:)=double(subs(dFs,[s t],[xBC(i) yBC(i)]));
end
T=null(FBC);  %-----(6)
Ks=T'*Kc*T  %-----(7)
Ms=T'*Mc*T  %-----(8)

% Eigenvalue solution
[PHI,EG]=eig(Ks,Ms);  %-----(9)
EgN=sort(diag(EG));
Wn=sqrt(EgN)

eigError=eigErr(Wn,lamDATA.ilam,bcType,a_n,b_n,abd,ipoltype,pDeg)



%function [ EgN,Ks,Ms,PHI,EG ] = PlateVibSymIso2(panXY,lamDATA,ipoltype,pDeg,plot_flag,isoltype,intMethod,bcMethod,bcType)
function [ EgN,Ks,Ms,PHI,EG ]=PlateVibSymIso2(panXY,lamDATA,solDATA)
%PlateVibSymIso: Solves plate vibration using symbolic notation -
%                in isoparametric coordinates
%   Detailed explanation goes here

display('============== Solve in {s,t} & {x,y} ===============');

[ipoltype,pDeg,Mg,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

% Symbolic variables
syms x y
syms s t
syms a b th       % plate length, width, thickness
syms rho nu D   % plate properties (D=Et^3/(1-nu^2))

% Symbolic polynomial basis function
[F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_Sym(ipoltype,pDeg) %-----(1)
F=subs(F, [x y], [s t]);
dFs=subs(dFx, [x y], [s t]);
dFt=subs(dFy, [x y], [s t]);
dFss=subs(dFxx, [x y], [s t]);
dFtt=subs(dFyy, [x y], [s t]);
dFst=subs(dFxy, [x y], [s t]);

% Gauss Legendre points and weights
if(intMethod == 2)
    order_1d=[pDeg+1,pDeg+1];
    order_nd = prod ( order_1d(1:2) );
    [gl_pts,gl_wts]=gl_grid( 2, order_1d, order_nd );
    if(plot_flag ==3)
        gl_grid_display2( order_1d, gl_pts );
    end
    sI=gl_pts(1,:);
    tI=gl_pts(2,:);
end

% Compute [Mc] and [Kc] in (s,t) coords
Kc=0;Mc=0;
for i=1:length(sI);
    psi=sI(i); eta=tI(i); gl_wt=gl_wts(i);
    [xI,yI]=ISO_st_to_xy(psi,eta,panXY);
    %[F,dFs,dFt,dFst,dFss,dFtt]=FdF_Sym(pDeg,'st');
    
    H=subs(F,[s,t],[xI,yI]); Fs=subs(dFs,[s,t],[xI,yI]); Ft=subs(dFt,[s,t],[xI,yI]);
    Fst=subs(dFst,[s,t],[xI,yI]); Fss=subs(dFss,[s,t],[xI,yI]); Ftt=subs(dFtt,[s,t],[xI,yI]);
    FssFtt_sym = Fss'*Ftt + Ftt'*Fss;

    J = jacob2D_Iso( psi, eta, panXY ); %??

%     Kc = Kc + D*(Fss.'*Fss + Ftt.'*Ftt + 2*nu*FssFtt_sym + 2*(1-nu)*Fst.'*Fst)*det(J)*gl_wt;
%     Mc = Mc + H.'*rho*th*H*det(J)*gl_wt;
    Kc = Kc + D*(Fss.'*Fss + Ftt.'*Ftt + 2*nu*FssFtt_sym + 2*(1-nu)*Fst.'*Fst)*det(J)*gl_wt;
    Mc = Mc + H.'*rho*th*H*det(J)*gl_wt;
end

% Get numerical values
[a_n, b_n, t_n, rho_n, D_n, nu_n, abd]=get_numeric_vals(lamDATA,panXY);
Mc=double(subs(Mc, [rho,th,a,b], [rho_n,t_n,a_n,b_n]))  %-----(4)
Kc=double(subs(Kc, [D,nu,a,b], [D_n,nu_n,a_n,b_n]))  %-----(5)

% BCs - clamped at x=0
numBC=pDeg+1;
yBC=linspace(0,b_n,numBC);
xBC=linspace(0,0,numBC);
IBC=0
for i=1:numBC
    IBC=IBC+1;
    FBC(IBC,:)=double(subs(F,[s t],[xBC(i) yBC(i)]));
    IBC=IBC+1;
    FBC(IBC,:)=double(subs(dFs,[s t],[xBC(i) yBC(i)]));
end
T=null(FBC);  %-----(6)
Ks=T'*Kc*T  %-----(7)
Ms=T'*Mc*T  %-----(8)

% Eigenvalue solution
[PHI,EG]=eig(Ks,Ms);  %-----(9)
EgN=sort(diag(EG));
Wn=sqrt(EgN)

%eigError=eigErr(Wn,lamDATA.ilam,bcType,a_n,b_n,abd,ipoltype,pDeg)
eigError=eigErr(panXY,solDATA,lamDATA,Wn);

%%
function [a_n, b_n, t_n, rho_n, D_n, nu_n, ABD]=get_numeric_vals(lamDATA,panXY)
% Panel geometry
a_n=panXY(2,1) - panXY(1,1);
b_n=panXY(3,2) - panXY(2,2);
%% Read laminate stack
%[imat, theta0, theta1, thk, mat_names] = readLam(ilam);

t_n=sum(lamDATA.thk);
%% Material properties
[e1,e2,nu12,g12,Rho,mtype] = get_mat_props(lamDATA.mat_names);
rho_n=Rho(1,1)/(12*32.2); nu_n=nu12(1,1);
%% Laminate [ABD]
%ABD=get_abd(a_n, b_n, 0, theta0, theta1 ,thk ,e1 ,e2 ,nu12 , g12);
%ABD=get_abd(a_n, b_n, 0, lamDATA);
ABD=get_abd(panXY, 0, 0, lamDATA);
D_n=ABD(4,4);



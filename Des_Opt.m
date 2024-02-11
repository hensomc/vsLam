function [XOPT,his]= Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% Design Optimization - Main driver function
disp('DES_OPT called');

plot_hist=1;
numX = length(desDATA.X0)

%desResults=cell(length(desDATA.Xstart));
desResults=[];

% Loop over start points
disp('length(desDATA.Xstart='); length(desDATA.Xstart)
for i=1:length(desDATA.Xstart)
    desDATA.X0 = desDATA.Xstart{i}

    switch desDATA.drespType
        case {'WEIGHT'}
            [XOPT,his]= Weight_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA);
        case {'EIGN'}
            [XOPT,his]= Vib_Freq_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA);
        case {'BUCKL'}
            [XOPT,his]= Buckl_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA);
        case {'STRAIN'}
            [XOPT,his]= Strain_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA);
        case {'WEIGHT_EIGN'}
            [XOPT,his,hisF]= Weight_Eign_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA);
            numX=length(XOPT);
            [his hisF]
            his=[his hisF(:,numX+1)]
        otherwise
            return;
    end
        
    %desResults{i} = {XOPT his(numX+1,end) his(numX+2,end)};
    %desResults = [desResults; XOPT(1:end) his(end,numX+1) his(end,numX+2)];
    desResults = [desResults; XOPT(1:end) his(end,numX+1:end)];

    % Plot results
    plotDes(XOPT,his,panXY,lamDATA,solDATA,ldsDATA,desDATA,plot_hist,i);
    
end
% celldisp(desResults)
% t1=cell2table(desResults);
% disp(t1)
format long
desResults

% Find optimum design values
I=find(max(max(desResults(:,numX+1))));
XOPT=desResults(I,1:numX);

disp('Pausing to allow visual check, hit any key to continue...')
pause


% Check optimal solution (Ritz and FEM)
solDATA.pDeg=8;
solDATA.Mg=solDATA.pDeg;
solDATA.femSol=1;
solDATA.plot_flag=1;


% Update laminate design values
if(lamDATA.tDeg>0)
    lamDATA.thk_coeff = XOPT;
    % Draw surface thickness constraint
    [g,heq]=vlam_Thk_Con(XOPT,panXY,lamDATA,solDATA);
else
    switch desDATA.drespType
        case {'WEIGHT' 'EIGN'}
            % lamDATA.thk = XOPT;  % constant thk
            lamDATA.thk = XOPT(1:2);
            
            %     lamDATA.theta0(1)=XOPT(3);  % fiber angles
            %     lamDATA.theta1(1)=XOPT(4);
            %     lamDATA.theta0(2)=-XOPT(3);
            %     lamDATA.theta1(2)=-XOPT(4);
            
        case {'BUCKL' 'STRAIN'}
            % Fiber angles only
            %setVarLamTheta( lamDATA, round(XOPT(1)), round(XOPT(2)) );
            
            nply=length(lamDATA.theta0);
            for i=1:nply/2
                if(mod(i,2) == 1)
                    lamDATA.theta0(i)=round(XOPT(1));
                    lamDATA.theta1(i)=round(XOPT(2));
                else
                    lamDATA.theta0(i)=-round(XOPT(1));
                    lamDATA.theta1(i)=-round(XOPT(2));
                end
            end
            j=1;
            for i=nply:-1:nply/2+1
                lamDATA.theta0(i)=lamDATA.theta0(j);
                lamDATA.theta1(i)=lamDATA.theta1(j);
                j=j+1;
            end
            
        otherwise
            return;
    end
end

% Compare Ritz optimized solution with Nastran model
[C,ritzEig,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
disp('finished vLam call');
[nasEig,nasDisp,nasEps,XYData1]=fe_mesh2D(panXY,solDATA,lamDATA,ldsDATA,solDATA.Mg+1,solDATA.Mg+1,solDATA.bcType);
disp('finished fe_mesh2D call');
getSolErr(panXY,solDATA,lamDATA,C,ritzEig,ritzDisp,ritzEps,sMesh,tMesh,nasEig,nasDisp,nasEps,XYData1);
disp('finished getSolErr call');

        
function [XOPT,his]= Weight_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% Weight_Des_Opt.m : Min weight design optimization w/ buckling constraint
disp('Min weight design optimization w/ buckling constraint');

global his;
his=[];

global hisG;
hisG=[];

% Handles to objective and constraint functions
OBJ=@Weight_Des_Obj;
CON=@gen_Con;

% Thickness inequality constraints
A=[];b=[];
if( lamDATA.tDeg > 0 )
    [A]=getInequalLayerCon(panXY,lamDATA,solDATA)
     tmax(size(A))= 0.150;
     tmin(size(A))= 0.001;

%   Min thk constraint
    A=-A';
    b=-tmin;
    
%   Max & Min thk constraint
%     A=[A';-A'];
%     b=[tmax;-tmin];
end

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
disp('XOPT=');
OP=optimset('disp','iter')  % Display optimization/iteration details
%XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA)
XOPT=fmincon(OBJ,X0,A,b,[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA,desDATA)

% Place design and constraint history into one matrix
his= [his hisG];



function [XOPT,hisF]= Vib_Freq_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_OPT.m : Max frequency design optimization w weight/vol constraint
disp('Optimize plate for max frequency subject to volume constraint');

global his
his=[];

% Handles to objective and constrint functions
OBJ=@vlam_Freq_Obj;
CON=@vlam_Volume_Con;

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0
%his = X0;  % Initial values X0 not captured by "his" - this fails out

% Inequality layer-thickness constraints
if( lamDATA.tDeg > 0 )
    [A]=getInequalLayerCon(panXY,lamDATA,solDATA)
    %pause
%     A = [A;-A];
     tmax(size(A))= 0.150;
     tmin(size(A))= 0.001;
%     b = [tmax;-tmin];

%   Min thk constraint
    A=-A';
    b=-tmin;
    
%   Max & Min thk constraint
%     A=[A';-A'];
%     b=[tmax;-tmin];
end

% Optimize
OP=[];
OP=optimset('disp','iter')  % Display optimization/iteration details
disp('XOPT=');
XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA)
%XOPT=fmincon(OBJ,X0,A,b,[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA)

function [XOPT,his]= Buckl_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_OPT.m : Max buckling eigenvalue design optimization w weight/vol constraint
disp('Optimize plate for max buckling eigenvalue subject to volume constraint');

global his
his=[];

% Handles to objective and constrint functions
OBJ=@Buckl_Obj;
CON=@fibK_Con;

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
OP=optimset('disp','iter')  % Display optimization/iteration details
disp('XOPT=');
XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA,desDATA)


function [XOPT,his]= Strain_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_OPT.m : Minimize Max strain design optimization w weight/vol constraint
disp('Minimize Max strain design optimization w weight/vol constraint');

global his
his=[];

% Handles to objective and constrint functions
OBJ=@Strain_Obj;
%CON=@fibK_Con;
CON=@vlam_Volume_Con;

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
OP=optimset('disp','iter')  % Display optimization/iteration details
disp('XOPT=');
XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA,desDATA)


%function [f,hx,xC]=Weight_Des_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
function [f]=Weight_Des_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% Weight_Des_Obj.m - weight design objective function g

global his;

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness
    lamDATA.thk=XDES(1:2); % constant thickness
%     lamDATA.theta0(1)=XDES(3);  % fiber angles
%     lamDATA.theta1(1)=XDES(4);
%     lamDATA.theta0(2)=-XDES(3);
%     lamDATA.theta1(2)=-XDES(4);
else
    lamDATA.thk_coeff=XDES; % variable thickness
end
f = panelWeightVthk(solDATA,lamDATA,panXY );

% Get current constraint value
%[g,heq]=buckl_Con(XDES,panXY,lamDATA,solDATA,ldsDATA);
% [C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
% g=-EgN(1)+1;

% Add XDES, Weight and eigenvalue to design history
%his=[his;XDES f g];
his=[his;XDES f ];


function [g,heq]=buckl_Con(XDES,panXY,lamDATA,solDATA,ldsDATA)
%% buckl_Con.m : buckling constraint function

global his

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness
    lamDATA.thk=XDES(1:2); % constant thickness
%     lamDATA.theta0(1)=XDES(3);  % fiber angles
%     lamDATA.theta1(1)=XDES(4);
%     lamDATA.theta0(2)=-XDES(3);
%     lamDATA.theta1(2)=-XDES(4);
else
    lamDATA.thk_coeff=XDES; % variable thickness
end

% Buckling eigenvalue solution
solDATA.isoltype=3;
[C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
g=-EgN(1)+1;

heq=0;


function [g,heq]=gen_Con(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% gen_Con.m : general constraint function

global hisG

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness
    lamDATA.thk=XDES(1:2); % constant thickness
%     lamDATA.theta0(1)=XDES(3);  % fiber angles
%     lamDATA.theta1(1)=XDES(4);
%     lamDATA.theta0(2)=-XDES(3);
%     lamDATA.theta1(2)=-XDES(4);
else
    lamDATA.thk_coeff=XDES; % variable thickness
end

% Static solution
g_strain=[];
solDATA.isoltype=3;
[C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
[g_strain,heq_strain]=strain_Con(ritzEps,panXY,lamDATA,solDATA,ldsDATA);

% Buckling eigenvalue solution
g_buckl=[];
% solDATA.isoltype=2;
% [C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
% g_buckl=-EgN(1)+1;

% Thickness constraints
g_thk=[];heq_thk=[];
if(lamDATA.tDeg>0)
    [g_thk,heq_thk]=vlam_Thk_Con(XDES,panXY,lamDATA,solDATA);
end

%g=[g_buckl; g_strain; g_thk];
g=[g_buckl; g_strain];
heq=0;

% Add constraint to history
%his_g_buckl=[his_g_buckl;g_bukl];
hisG=[hisG; max(max(g))];


%function [f,hx,xC]=vlam_Freq_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
function [f]=vlam_Freq_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_Obj.m : plate frequency objective function

global his
global hisF

% original
%lamDATA.thk=XDES;
%lamDATA.thk_coeff=[1 0 0 0];

% lamDATA.thk=XDES(1:2);
% lamDATA.thk_coeff=[XDES(3) 0 0 0];
% lamDATA.thk_coeff=[XDES(1) 0 XDES(2) 0]; % t=linear with x
% lamDATA.thk_coeff=[XDES(1) XDES(3) XDES(2) XDES(4)]; % t=bi-linear with x & y
if(lamDATA.tDeg==0)
    lamDATA.thk=XDES; % constant thickness
else
    lamDATA.thk_coeff=XDES; % variable thickness
end


% Eigenvalue solution
solDATA.isoltype=1;
[C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);

% Set objective function value: maximize 1st frequency
WN=sqrt(EgN);
f=-WN(1);

% Add current weight and thicknesses to design history his
%Wt = panelWeight(lamDATA.thk,lamDATA.rho_ply,panXY);
Wt = panelWeightVthk(solDATA,lamDATA,panXY );
%his=[his;XDES WN(1) Wt];
% his=[his;XDES(1:2) WN(1) Wt];
%his=[his;XDES WN(1)];
hisF=[hisF;XDES WN(1)];


function [f,hx,xC]=Buckl_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% Buckl_Obj.m : plate buckling eignevalue objective function

global his

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness

    % Fiber angles only    
    %setVarLamTheta( lamDATA, XDES(1), XDES(2) );    
    
    nply=length(lamDATA.theta0);
    for i=1:nply/2
        if(mod(i,2) == 1)
            lamDATA.theta0(i)=XDES(1);
            lamDATA.theta1(i)=XDES(2);
        else
            lamDATA.theta0(i)=-XDES(1);
            lamDATA.theta1(i)=-XDES(2);
        end
    end
    j=1;
    for i=nply:-1:nply/2+1
        lamDATA.theta0(i)=lamDATA.theta0(j);
        lamDATA.theta1(i)=lamDATA.theta1(j);
        j=j+1;
    end
else
    %lamDATA.thk_coeff=XDES; % variable thickness
end

% Eigenvalue solution
[C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);

% Set objective function value
f=-EgN(1);

% Add constraint to design history his
Wt = panelWeightVthk(solDATA,lamDATA,panXY );
[g,heq]=fibK_Con(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);

%his=[his;XDES EgN(1) Wt];
his=[his;XDES EgN(1) g];


function [f,hx,xC]=Strain_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% Strain_Obj.m : strain objective function

global his

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness

    % Fiber angles only    
    %setVarLamTheta( lamDATA, XDES(1), XDES(2) );    
    
    nply=length(lamDATA.theta0);
    for i=1:nply/2
        if(mod(i,2) == 1)
            lamDATA.theta0(i)=XDES(1);
            lamDATA.theta1(i)=XDES(2);
        else
            lamDATA.theta0(i)=-XDES(1);
            lamDATA.theta1(i)=-XDES(2);
        end
    end
    j=1;
    for i=nply:-1:nply/2+1
        lamDATA.theta0(i)=lamDATA.theta0(j);
        lamDATA.theta1(i)=lamDATA.theta1(j);
        j=j+1;
    end
else
    %lamDATA.thk_coeff=XDES; % variable thickness
end


% Ritz solution
[C,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);

% Set objective function value: minimize max strain
ritzMax = [max(max(ritzEps(:,1))); max(max(ritzEps(:,2))); max(max(ritzEps(:,3))); max(max(ritzEps(:,4))); max(max(ritzEps(:,5)))];
ritzMin = [min(min(ritzEps(:,1))); min(min(ritzEps(:,2))); min(min(ritzEps(:,3))); min(min(ritzEps(:,4))); min(min(ritzEps(:,5)))];

% Get max strain
ritzMaxAbsAllComp = max(max(abs(ritzMax)));
ritzMinAbsAllComp = max(max(abs(ritzMin)));
ritzMaxAbs = max([ritzMaxAbsAllComp ritzMinAbsAllComp]);

f=ritzMaxAbs;

% Add constraint to design history his
Wt = panelWeightVthk(solDATA,lamDATA,panXY );
[g,heq]=fibK_Con(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);

%his=[his;XDES EgN(1) Wt];
his=[his;XDES f g];

function [g,heq]=strain_Con(ritzEps,panXY,lamDATA,solDATA,ldsDATA)
%% Strain_Obj.m : strain objective function

% Allowables
allow_epsx = .0018;
allow_epsy = .0018;
allow_epsxy= .0027;
allow_epsxz= .0015;
allow_epsyz= .0015;

% Find Min/Max strains
ritzMax = [max(max(ritzEps(:,1))); max(max(ritzEps(:,2))); max(max(ritzEps(:,3))); max(max(ritzEps(:,4))); max(max(ritzEps(:,5)))];
ritzMin = [min(min(ritzEps(:,1))); min(min(ritzEps(:,2))); min(min(ritzEps(:,3))); min(min(ritzEps(:,4))); min(min(ritzEps(:,5)))];

% Get max strain
ritzMaxAbsAllComp = max(max(abs(ritzMax)));
ritzMinAbsAllComp = max(max(abs(ritzMin)));
ritzMaxAbs = max([ritzMaxAbsAllComp ritzMinAbsAllComp]);

% Calc constraint
g=ritzMaxAbs/0.001800 -1;

% Fiber Strains: Eps-x
g_epsx = ritzEps(:,1)/allow_epsx -1;

% Fiber Strains: Eps-y
g_epsy = ritzEps(:,2)/allow_epsy -1;

% Shear Strains: Eps-xy
g_epsxy = ritzEps(:,3)/allow_epsxy -1;

% Transverse Strain
g_epsxz = ritzEps(:,4)/allow_epsxz -1;

% Transverse Strain
g_epsyz = ritzEps(:,5)/allow_epsyz -1;

% Combine constraints
g=[g_epsx; g_epsy; g_epsxy; g_epsxz; g_epsyz];

heq=0;  % Eq constraints



function [g,heq]=vlam_Volume_Con(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Volume_Con.m - plate volume constraint function g

% Compute volume constraint value
% lamDATA.thk=XDES;
% lamDATA.thk=XDES(1:2);

if(lamDATA.tDeg==0)
    lamDATA.thk=XDES; % constant thickness
else
    lamDATA.thk_coeff=XDES; % variable thickness
end
%Wt = panelWeight(lamDATA.thk,lamDATA.rho_ply,panXY);
Wt = panelWeightVthk(solDATA,lamDATA,panXY );
WtBound = 2.0; % 1.0 lbs max weight
g=Wt/WtBound-1; %Ineq constraints

if(lamDATA.tDeg > 0)
    % Thickness constraints
    %[g_thk,heq_thk]=vlam_Thk_Con(XDES,panXY,lamDATA,solDATA);
    
    % Combine panel weight constraint and thickness constraints into vector
    %g = [g; g_thk];
end

heq=0;  % Eq constraints



function [g,heq]=fibK_Con(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% fibK_Con.m - fiber curvature constraint function

% Allowable curvature
Kallow=1/desDATA.Rallow;

% Start by checking only 1st ply
T0=XDES(1);
T1=XDES(2);

a=abs(panXY(2,1)-panXY(1,1));
x=linspace(0,a/2);

fibK = getFibCurvature(a,T0,T1,x);
fibK_min=min(fibK);

% Find curvature less than input allowable
I=find(fibK<=Kallow);

% Normalize constraint
%g=fibK/Kallow-1
%g=fibK_min/Kallow-1;
g=Kallow/fibK_min-1;

heq=0;  % Eq constraints


function [g,heq]=vlam_Thk_Con(XDES,panXY,lamDATA,solDATA)
%% vlam_Thk_Con.m - plate thickness constrain function g

% Enforces postive thickness over plate surface

% Set lamDATA thk coeffs to optimized values
lamDATA.thk_coeff=XDES;

if(lamDATA.tDeg > 0)
    % mesh data for thk evaluation pts
    [ meshDATA, intDATA ] = meshGen( panXY, solDATA );
    [psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts] = getIntData(intDATA);

    % compute surface representation of thickness
%     thk_coeff_mat = repmat(lamDATA.thk_coeff,20,1);
%     z1 = thk_legendre(lamDATA.tDeg,psiIMesh,etaIMesh)*thk_coeff_mat';
%     disp('size of lamDATA.thk_coeff='); size(lamDATA.thk_coeff)
%     disp('size of thk_coeff_mat='); size(thk_coeff_mat)
%     disp('size of z1='); size(z1)
%     disp('size of psiIMesh='); size(psiIMesh)
%     disp('size of etaIMesh='); size(etaIMesh)
%    surface(psiIMesh,etaIMesh,z1,'FaceColor',[0.5 1.0 0.5], 'EdgeColor', 'none');

    dim = length(psiIMesh)^2;
    psiIMesh2=reshape(psiIMesh,1,dim); etaIMesh2=reshape(etaIMesh,1,dim);
    z1 = thk_legendre(lamDATA.tDeg,psiIMesh2,etaIMesh2)*lamDATA.thk_coeff';
    z1 = reshape(z1,size(psiIMesh));
    %surface(psiIMesh,etaIMesh,z1,'FaceColor',[0.5 1.0 0.5], 'EdgeColor', 'none');
    
    % define a plane at z=0
    z2 = zeros(size(psiIMesh));
    %surface(psiIMesh,etaIMesh,z2,'FaceColor',[1.0 0.5 0.0], 'EdgeColor', 'none');

    % intersection of computed thk surface with min-thk surface
    zdiff = z1 - z2;
    zdiff_norm=zdiff/0.125 - 1.0; % arbitrary set Zupper-Zlower=0.125
    %test_neg=all(zdiff < 0)
    
    % return values
    g=reshape(zdiff_norm,length(zdiff)^2,1); 
    heq=0;
    
    % get min-max thicknesses
    tmax = max(max(zdiff));
    tmin = min(min(zdiff));

    % Draw thickness surface and contour of intersection
    if(solDATA.plot_flag > 0)
        figure
        surface(psiIMesh,etaIMesh,z1,'FaceColor',[0.5 1.0 0.5], 'EdgeColor', 'none');
        surface(psiIMesh,etaIMesh,z2,'FaceColor',[1.0 0.5 0.0], 'EdgeColor', 'none');
        
        %  Contours of surface intersection
        C = contours(psiIMesh, etaIMesh, zdiff, [0 0]);
        %disp('Size of contour intersections = '); C
        %disp('Norm of C = ');norm(C)
        
        % Extract the x- and y-locations from the contour matrix C.
        xL = C(1, 2:end);
        yL = C(2, 2:end);
        
        % Interpolate on the first surface to find z-locs for the intersection
        zL = interp2(psiIMesh, etaIMesh, z1, xL, yL);
        
        % Visualize the line.
        line(xL, yL, zL, 'Color', 'k', 'LineWidth', 3);     
    end
end

function [A]=getInequalLayerCon(panXY,lamDATA,solDATA)
%% vlam_Thk_Con.m - plate thickness constrain function g

% Inequality constraint Matrix A - enforces postive thickness over plate surface
if(lamDATA.tDeg > 0)
    % mesh data for thk evaluation pts
    [ meshDATA, intDATA ] = meshGen( panXY, solDATA );
    [psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts] = getIntData(intDATA);

%     dim = length(psiIMesh)^2;
%     psiIMesh2=reshape(psiIMesh,1,dim); etaIMesh2=reshape(etaIMesh,1,dim);
%     A = thk_legendre(lamDATA.tDeg,psiIMesh2,etaIMesh2)';

    [psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(20);
    dim = length(psiBCMesh)^2;
    psiBCMesh2=reshape(psiBCMesh,1,dim); etaBCMesh2=reshape(etaBCMesh,1,dim);
    A = thk_legendre(lamDATA.tDeg,psiBCMesh2,etaBCMesh2)';

end



function [XOPT,his,hisF]= Weight_Eign_Des_Opt(panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% Weight_Des_Opt.m : Min weight design optimization w/ strain constraint
disp('Min weight design & max freq optimization w/ strain constraint');

global his
his=[];

global hisF
hisF=[];

global hisG
hisG=[];

% Handles to objective and constraint functions
OBJ=@vib_weight_obj; % must return a vector
CON=@gen_Con;

wfactor=[1;1];
goal=[1.0;100];

% Thickness inequality constraints
A=[];b=[];
if( lamDATA.tDeg > 0 )
    [A]=getInequalLayerCon(panXY,lamDATA,solDATA)
     tmax(size(A))= 0.150;
     tmin(size(A))= 0.001;

%   Min thk constraint
    A=-A';
    b=-tmin;
    
%   Max & Min thk constraint
%     A=[A';-A'];
%     b=[tmax;-tmin];
end

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
disp('XOPT=');
OP=optimset('disp','iter')  % Display optimization/iteration details
%XOPT=fgoalattain(OBJ,X0,goal,wfactor,A,b,[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA,desDATA)
XOPT=fgoalattain(OBJ,X0,goal,wfactor,A,b,[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA,desDATA)

% Place design and constraint history into one matrix
%his= [his hisG];

%function [weight_obj, vib_freq_obj] = vib_weight_obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
function [out] = vib_weight_obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
%% vibration and weight objective function
weight_obj=Weight_Des_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);
vib_freq_obj=vlam_Freq_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);

out=[weight_obj; vib_freq_obj];





function plotDes=plotDes(XOPT,his,panXY,lamDATA,solDATA,ldsDATA,desDATA,plot_hist,istart)
%% Plot design solutions

disp('Plotting design history');

% his = history matrix (design variables, obj fn, panel weight)
%       his(1:2) = design variables x(1), x(2)
%       his(3)   = obj fn
%       his(4)   = panel weight

format long;
HIS=his

numX=length(XOPT);

if(plot_hist > 0)
    figure
    % Design variable history
    subplot(3,1,1)
    %plot(his(:,1:2),'linewidth',3)
    %plot(his(:,1:4),'linewidth',3)
    plot(his(:,1:numX),'linewidth',3)
    title('\bfDesign history')
    ylabel('\bfDesign variables')
    
    % Objective funtion history
    subplot(3,1,2)
    %plot(his(:,3),'linewidth',3)
    %plot(his(:,5),'linewidth',3)
    plot(his(:,numX+1),'linewidth',3)
    hold on;
    if(strcmp(desDATA.drespType,'WEIGHT_EIGN'))
        plot(his(:,numX+2),'linewidth',3)
        hold off;
        
        % Pareto Front
        subplot(3,1,3)
        plot( his(:,numX+1), his(:,numX+2),'-o','linewidth',2 )
        
        if(lamDATA.tDeg > 0)
            solDATA.plot_flag=1;
            [g,heq]=vlam_Thk_Con(XOPT,panXY,lamDATA,solDATA);
        end

        return;
    end
    ylabel('\bfObjective')
    
    % Constraint function history
    subplot(3,1,3)
    %plot(his(:,4),'linewidth',3)
    %plot(his(:,6),'linewidth',3)
    plot(his(:,numX+2),'linewidth',3)
    ylabel('\bfConstraint')
    xlabel('\bfFunction evaluation number')

    print('design_history', '-dpng');

end

% Compute data for objective function contours
global his
his=[];

global hisG
hisG=[];

% Range of design variables for Obj Fn contours
x1=linspace(desDATA.XL(1),desDATA.XU(1)+0.1,10);
x2=linspace(desDATA.XL(2),desDATA.XU(2)+0.1,10);

% Early attempt to use polynomial coeffs
% x1=linspace(-1.0,1.0,10);
% x2=linspace(-1.0,1.0,10);
% x3=linspace(-1.0,1.0,10);
% x4=linspace(-1.0,1.0,10);


%% Contour plot of objective function and constraints
if(lamDATA.tDeg==0)
    if(istart ==1)  % only do this once if mult. start pts
        [xx,yy]=meshgrid(x1,x2);
        for i=1:length(x1)
            for j=1:length(x2)
                XDES=[x1(i) x2(j)]; %original
                
                %XDES=[x1(i)  x2(j) thk_coeff(1)];
                %XDES=[x1(i)  x2(j)]; % linear with x
                %XDES=[x1(i)  x2(j)  x3(j)  x4(j)];  % bi-linear with x,y
                %XDES=[x1(i)  0.00  0.00  x2(j)  0.00  0.00   x3(j)  0.00  0.00];  % bi-linear with x,y
                
                % compute obj fn value
                switch desDATA.drespType
                    case {'WEIGHT'}
                        f=Weight_Des_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);
                        g=gen_Con(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);
                    case {'EIGN'}
                        f=vlam_Freq_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA);
                    case {'BUCKL'}
                        f=-Buckl_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);
                    case {'STRAIN'}
                        f=Strain_Obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA);
                    otherwise
                        f=0;
                        return;
                end
            end
        end
        his = [his hisG];
        CHIS=his;  % Note: 'his' is updated automatically from Obj Fn calls above
        
        % Objective Fn Contours: OBJC
        figure
        OBJC=reshape(his(:,numX+1),size(xx));  %his(:,numX+1) = obj fn
        objmin=min(min(OBJC));objmax=max(max(OBJC));
        [cc,hh]=contour(xx,yy,OBJC);  % transposed OBJC
        set(hh,'linewidth',3)
        clabel(cc)
        hold on
        
        %     % ========================================================
        %     % Draw feasible design region based on min steering radius
        %     % ========================================================
        %     a=abs(panXY(2,1)-panXY(1,1));
        %     x=linspace(0,a/2);
        %     xx1=linspace(desDATA.XL(1),desDATA.XU(1),25);
        %     xx2=linspace(desDATA.XL(2),desDATA.XU(2),25);
        %     [xxc,yyc]=meshgrid(xx1,xx2);
        %
        %     for i=1:length(xx1)
        %         for j=1:length(xx2)
        %             fibK = getFibCurvature(a,xx1(i),xx2(j),x);
        %             fibK_min(i,j) = min(min(abs(fibK)));
        %             fibK_max(i,j) = max(max(abs(fibK)));
        %         end
        %     end
        %     % Find curvatures less than input allowable
        %     I=find(fibK_min<=1/desDATA.Rallow);
        %
        %     % Plot allowable fiber curvature points
        %     X=xxc(:);Y=yyc(:);
        %     plot(X(I),Y(I),'y.','MarkerSize',50);
        %     % ========================================================
        
        % Constraint Contours: CONC
        CONC=reshape(his(:,numX+2),size(xx))  %his(:,numX+2) = constraint
        gmax=max(max(CONC))
        gmin=min(min(CONC))
        if(gmax ~= gmin)
            if(gmax > 10)
                gmax=1;
            end
            CL=linspace(gmin,gmax,10);
            CL = [CL 0]; % ensure that pt for feasible constraint boundary included in level
            CL = sort(CL);
            %CL=linspace(gmin,gmax,5);
            [cc1,hh1]=contour(xx,yy,CONC,CL);clabel(cc1)    % all contours
            
            %[cc1,hh1]=contour(xx,yy,CONC,[0 0]);clabel(cc1)  % constraint boundary = [0 0]
            set(hh1,'linewidth',2,'color','k')
            
        else
        end
        
        %     % Redraw obj contours to overlay on feasible region
        %     [cc,hh]=contour(xx,yy,OBJC);  % transposed OBJC
        %     set(hh,'linewidth',3)
        %     clabel(cc)
        
    else
        hold on
    end
    
    % Overlay optimized result values
    plot(XOPT(1),XOPT(2),'rp','linewidth',5)
    plot(HIS(:,1),HIS(:,2),'b--','linewidth',3)
    
    % plot(XOPT(1),XOPT(4),'rp','linewidth',5)
    % plot(HIS(:,1),HIS(:,4),'b--','linewidth',3)
    
    % title({'\bfObjective Function and Constraint Contours'; ...
    %        '\bfTwo-Layer Cantilevered Plate'; ...
    %        '\bfAluminum/CarbonEpoxy'});
    [panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);
    
    title({'\bfObjective Function and Constraint Contours'; ...
        panTitle; lamTitle});
    xlabel('\bfDesign Variable X_1');
    ylabel('\bfDesign Variable X_2');
    
    hold off
    print('design_space_coutours', '-dpng');
    
else
    %Draw panel thickness function
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(solDATA.pDeg,solDATA.intMethod,panXY,solDATA.plot_flag);
    lamDATA.thk_coeff=XOPT;
    layer_thk = getLayerThk(lamDATA.tDeg,psiI,etaI,lamDATA,1);
end

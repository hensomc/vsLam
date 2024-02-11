%% function  [XYData,EData]= M_Mesh_2D_B_May7_2015(XYData,Region,ElemType,IPLOT) % Generalization of  M_Mesh
%
%% Purpose:
% The function M_Mesh2D_A can be used to generate meshes for 2-D FEA
% for a composite domain
%% Usage:
%  [XYData,EData]= M_Mesh_2D_A(XYData,Region,ElemType,IPLOT)
%% Input: 
%% 1. Key-points coordinates : these are key points of the domain
%      XYData =[ Key_Point_No   X    Y];
%% 2. Define regions use key points, each region per line
% Region=[  Region_NO   K1 K2  K3  K4  N12   N14];
%     Region_No= Region Number of a quadrilateral region to be meshed
%     I1 to I4 = key points ID for the region, must be
%                counterclockwise
%     N12 = number of nodes to be defines from key points K1 to K2
%     N14 = number of nodes to be defines from key points K1 to K4
%  
%% 3. Element type :
% ElemType='T3'; %  for 3-noded triangular element
% ElemType='Q4'; %  for 4-noded quadrilateral element
% ElemType='Q8'; %  for 8-noded quadrilateral element
%% Special case: one rectangular domain
% for the special case of a rectangular daomain,
% The model is defined by 4 keynodes and 1 region. 
% Then the mesh can be generated using 'T3', 'Q4', ]Q8' elements
% and 'R4','R8' and 'TU', 'TD' elements.
%% Written by B.P. Wang, Dec 6,2014
%
%(c) B.P. Wang, 2014,2015

function [XYData,EData]=M_Mesh2D_B_May7_2015(XYData,Region,ElemType,IPLOT)
%function  [XYData,EData]= M_Mesh_2D_A(XYData,Region,ElemType,IPLOT) % Generalization of  M_Mesh

%% Rectangular domain
x=XYData(:,2);
    XI=x(1);XF=x(2);
    NX=Region(1,6);
    xx=linspace(XI,XF,NX);
y=XYData(:,3);
    YI=y(1);YF=y(4);
    NY=Region(1,7);
    yy=linspace(YI,YF,NY);


if upper(ElemType(1:2))=='R4';
    [p,t,Elec,NodeXY]=RECT_4_MeshN(xx,yy,IPLOT);
NN=(1:length(p(:,1)))';
XYData=NodeXY;

end
if upper(ElemType(1:2))=='R8';
    [p,t,Elec,NodeXY]=RECT_8_MeshN(xx,yy,IPLOT);
XYData=NodeXY;
end


if upper(ElemType(1:2))=='TU';
    [XYData,Elec,XX,YY]=Tri_Mesh1_RECTN(xx,yy,IPLOT);end
if upper(ElemType(1:2))=='TD';
    [XYData,Elec,XX,YY]=Tri_Mesh2_RECTN(xx,yy,IPLOT);end
% if upper(ElemType(1:2))=='TX'; do not work, 12-8-14
%     [XYData,Elec,XX,YY]=Tri_Mesh_Cross2N(xx,yy,IPLOT);end


% %% Circular ring domain  to be done
% 
% 
% angles=x;radius=y;
% if ElemType=='RB';
%     [p,t,Elec,X,Y]=RECT_MeshN_Ring(angles,radius,IPLOT);
%     XTData=[X Y];EData=Elec;
% end
% 
% % if ElemType=='RIT3';
% %     [XYData,EData]= M_Mesh_T3_Function(XYData,Region)
% % end
% 
% if ElemType=='RC';
%  %   [XYData,EData]= M_Mesh_Q8_Function(XYData,Region)
%    [XYData,EData]=RECT_8_MeshN_Ring(x,y,IPLOT);
% end
% 
% 




%% Composite Domain
%x=XYData;y=Region;
ElemType=upper(ElemType);
if ElemType=='Q4';
    [XYData,EData]= M_Mesh_Q4_Function(XYData,Region);
    Elec=EData;
end

if ElemType=='T3';
    [XYData,EData]= M_Mesh_T3_Function(XYData,Region);
    Elec=EData;
end
% disp('Line 59')
if ElemType=='Q8';
    [XYData,EData]= M_Mesh_Q8_Function(XYData,Region);
    Elec=EData;
end

EData=Elec;


%% --------------------  end


%%  sub-functions ================================= 

%ISO_st_to_xy.m : Given (s,t) compute (x,y)
%function [x,y]=ISO_st_to_xy(s,t,XY)
% Input:  
%    s = s-coordinate of evaluation point
%    t = t-coordinate of evaluation point
%  XY = nodal coordinates=[x y]
%OUTPUT;
%    x = x-coordinate of evaluation point
%    y = y-coordinate of evaluation point
%    N = shape function at this point

function [x,y]=ISO_st_to_xy(s,t,XY)

NN=length(XY(:,1));
    N=N_ISO48(s,t,NN);
    x=N*XY(:,1);
    y=N*XY(:,2);
   
   

function  [XYData,EData]= M_Mesh_Q4_Function(XYData,Region)


%% Plot model regions
ITEXT=2;
figure;SeeRecModel(XYData,Region,ITEXT)
ITEXT=1;hold on;SeeRecModel(XYData,Region,ITEXT)
title('\bfRegions for mesh generation')
X=XYData(:,2);Y=XYData(:,3);
%% Number of elements nd nodes per region
NRegion=prod(Region(:,6:7)'); % Number of nodes in region
NN=cumsum(NRegion);
ERegion=prod([(Region(:,6)-1) (Region(:,7)-1)]'); % Number of nodes in region
NEN=cumsum(ERegion);
                        %disp('Line 14 =================================')
TNodeXY=[]; % with duplicate coordinates
tALL=[];
%% Mesh all regions
    XID=[]; % region ID, 5-7-2015
for iR=1:length(Region(:,1));
Ns=Region(iR,6);S=linspace(-1,1,Ns);
Nt=Region(iR,7);T=linspace(-1,1,Nt);
%figure % Plot domain
IPLOT=0;
[p,t,Elec]=RECT_MeshN(S,T,IPLOT);
%disp('line 159'),size(p),size(t),
nt=length(t(:,1));
% t=[t ones(nt,1)*iR]
XID=[XID;ones(nt,1)*iR];
% Elec=[Elec ones(nt,1)*iR]
% pause
%  disp('Line 24 ================================='),pause

%% Convert (s,t) into (x,y)
JABCD=Region(iR,2:5);
XE=X(JABCD);
YE=Y(JABCD);
XY=[XE YE];
ss=p(:,1);tt=p(:,2);
NodeXY=[];
for is=1:length(ss)
[x,y]=ISO_st_to_xy(ss(is),tt(is),XY);
NodeXY=[NodeXY;x y];
end
II=1:NRegion(iR);
if iR>1
II=(1:NRegion(iR))+NN(iR-1);
end
NodeTemp=[II' NodeXY];
TNodeXY=[TNodeXY;NodeTemp]; % with duplicate coordinates
if iR==1;nadd=0;end
if iR>1;nadd=NN(iR-1);end
tALL=[tALL;t+nadd];         % with node number to change
end

%% Define ElecTemp
IE=(1:NEN(length(Region(:,1))))';
ElecTemp=[IE tALL];
%disp('Line 185'),pause


%% Get rid of duplicated nodes
%% M_Mesh_Development_2.m
%% Get rid of duplicated nodes
XYT=TNodeXY;
format compact
XYKeep=[XYT(1:NN(1),[2 3])]; % keep all nodes in Region I
ILast=NN(1);           % Last node No.
NodeNew=(1:NN(1))';
for in=NN(1)+1:max(NN)
%for in=19:21
    XX=XYKeep(:,1);
    YY=XYKeep(:,2);
nK=length(XX);
   
     Xn=XYT(in,2);
     Yn=XYT(in,3);
  
     DX=(abs(Xn-XX));
     DY=(abs(Yn-YY));
    
     NI=find(DX+DY<1e-6);

     if length(NI)==0;XYKeep=[XYKeep;Xn Yn];
         ILast=ILast+1; NodeNew(in)=ILast;end
        
     if length(NI)>0;NodeNew(in)=NI;
     end
 end

%% Filtered nodes
IN=(1:length(XYKeep(:,1)))';
NodeXYNew=[IN XYKeep];
%% Change element connectivity data

J=tALL(:,1);JA=NodeNew(J);
  J=tALL(:,2);JB=NodeNew(J);
   J=tALL(:,3);JC=NodeNew(J);
    J=tALL(:,4);JD=NodeNew(J);
    IE=(1:max(NEN))';
%     size(IE)
%     size(JA)
%     II
%     pause
%    ELECNew=[IE JA JB JC JD  ];  % II=region number  5-7-2015
    ELECNew=[IE JA JB JC JD  XID];  % II=region number  5-7-2015
%% Plot mesh
XYData=NodeXYNew;
EData=ELECNew;
% %% write node number:
% ITEXT=1;
% SeeRecModel(XYData,EData,ITEXT)

%% write element number:
ITEXT=2;
 if max(NN)>30;ITEXT=0;end

figure
SeeRecModel(XYData,EData,ITEXT)
%% --------------------  end



%diary RECT_8_MeshN_Ring.m : do not wiok
%type RECT_8_Mesh

%% function [p,t,Elec,NodeXY]=RECT_8_MeshN(x1,x2,nx,y1,y2,ny,IPLOT);
%% Purpose:
%   Generate uniform rectangular mesh for 8-noded rectangular elements in a rectangular domain
%% Input
% x1,x2,nx   = minimum x, maximun x, number of points
% y1, y2, ny = minimum y, maximun y, number of points
%   IPLOT =1; % plot medh
%         =0, do not plot
%% Output
%   p  = [x y] nodal coordinates
%   t  = [JA JB JC JD], element connectivity
%  Elec= [ElenentNo JA JB JC JD JE JF JG JH], element connectivity
%  NodeXY= [NodeNo X  Y], element connectivity
%
%% B.P. Wang, July 19,2012

%% function [p,t,Elec,NodeXY]=RECT_8_MeshN_Ring(x,y,IPLOT);
function [NodeXY,Elec]=RECT_8_MeshN_Ring(x,y,IPLOT);

% Developed in REC8_Mesh_Development.m
[p,t,Elec]=RECT_MeshN_Ring(x,y,IPLOT);

% Expand ranges
VV=axis;
RX=VV(2)-VV(1);RY=VV(4)-VV(3);
RR=max([RX RY]);
a=0.1;
V=VV+[-a a -a a]*RR; axis(V)


%aa=0.2;AxisExpand(aa)
%% ===================================================
%% Add vertical nodes
%% Unidorm
% DY=(y2-y1)/(ny-1);  %% to be done
%     x=linspace(x1,x2,nx);
%     y=linspace(y1+DY/2,y2-DY/2,ny-1);
%     [X,Y]=meshgrid(x,y);
% px=X(:);py=Y(:);
% pp=[p;px py]
%% Non-Unidorm
ny=length(y);
nx=length(x);
DY=y(2:ny)-y(1:(ny-1));  %% to be done

ya=y(1:(ny-1))+DY/2;
     [X,Y]=meshgrid(x,ya);
 px=X(:);py=Y(:);
 pp=[p;px py]



plot(px,py,'ro','markersize',12)
nv=nx*(ny-1);
np=nx*ny+(1:nv);
text(px,py,int2str(np'),'color','r','Fontsize',15)
%% Add horizonztal nodes
%% Uniform
% DX=(x2-x1)/(nx-1);
%     y=linspace(y1,y2,ny);
%     x=linspace(x1+DX/2,x2-DX/2,nx-1);
%     [X,Y]=meshgrid(x,y);
%% Non-Uniform

 DX=x(2:nx)-x(1:(nx-1));

 x=x(1:nx-1)+DX/2;
 [X,Y]=meshgrid(x,y);
px=X(:);py=Y(:);
pp=[pp;px py]
plot(px,py,'rs','markersize',12)
nh=ny*(nx-1);
np=nx*ny+nv+(1:nh);
text(px,py,int2str(np'),'color','r','Fontsize',15)

%% Nodes H and F
nn=nx*ny;
nnx=nx*(ny-1);
nny=ny*(nx-1);
HH=nn+(1:(nx-1)*(ny-1))';
FF=HH+(ny-1);
[HH FF]

%% Nodes E
EE=[];
for ix=1:nx-1
    e1=(nn+nnx)+(ix-1)*ny+(1:(ny-1))';
    EE=[EE;e1];
end
%% Nodes G
GG=EE+1;
t58=[EE FF GG HH];
t8=[t t58];
NE=length(t8(:,1));
%% Output
Elec=[(1:NE)' t8 ones(NE,2)];
II=(1:length(pp(:,1)))';
NodeXY=[II pp];




%diary  Tri_Mesh_Cross2N.m
%type Tri_Mesh_Cross2

% %diary aa.m
%% function [XYData,Elec,XX,YY]=Tri_Mesh_Cross2(x1,x2,NDx,y1,y2,NDy,IPLOT)
 function [XYData,Elec,XX,YY]=Tri_Mesh_Cross2N(x,y,IPLOT)
%% In :  F:\00-FE-Book-Work\00-FE-Book-BP-added\Daily-work\2D-Mesh-6-3-12

% type  rect_mesh.m
% 
% %RECT_Mesh.m : Generate rectangular mesh
% function [p,t,Elec,X,Y]=RECT_Mesh(x1,x2,nx,y1,y2,ny,IPLOT)
% %Input

% ny=2*NDy+1;
nx=length(x);ny=length(y);
NDx=(nx-1)/2;
NDy=(ny-1)/2;
        if NDx==0.5;NDx=1;end
        if NDy==0.5;NDy=1;end

%     x=linspace(x1,x2,nx);
%     y=linspace(y1,y2,ny);
    [X,Y]=meshgrid(x,y);XX=X;YY=Y;
%% Dertermine center node number
ic=0;
t=[];
for ix=1:NDx
    for iy=1:NDy
        ic=ic+1;
        
        Na=2*(ix)*ny-ny;
        NC(ix,iy)=Na+(iy-1)*2+2;
        CN=Na+(iy-1)*2+2;
        CC(ic)=CN;
        A1=CN-(1+ny);
        A2=A1+1;
        A3=A2+1;
        
        A4=CN-1;
        A6=CN+1;
         A8=CN+ny;
         A7=A8-1;A9=A8+1;
         
         EL=[ CN A2 A1
             CN A1 A4
             CN A4 A7
             CN A7 A8
             CN A8 A9
             CN  A9 A6
             CN  A6 A3
             CN A3 A2];
     t=[t;EL];    
    end
end

        
        
        
        
    
%     t=t,disp('Line 60'),pause
    
    
    
NE=length(t(:,1));
IE=(1:NE)';
Elec=[IE t];

p=[X(:) Y(:)];
t=Elec(:,2:4);
IN=(1:length(X(:)))';
XYData=[IN p];
NElem=length(Elec(:,1));
    figure,plot(X(:),Y(:),'ro')
    text(X(:),Y(:),int2str((1:prod(size(X)))'))
hold on

for ie=1:NElem
    tN=t(ie,:);
    xye=sum(p(tN,:));
    xy=sum(p(Elec(ie,2:4),:))/3;
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
        xye=p(tN,:);
 
        
%    if IPLOT>0
    plot(xye(:,1),xye(:,2),'b','linewidth',2)
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
  %end

end
ICC=CC;
plot(X(ICC),Y(ICC),'bs','linewidth',3)





%ISO_xy_to_st.m Given (x,y) compute (s,t)
%function [s,t]=ISO_xy_to_st(x,y,XY)
% Input:  
%    x = x-coordinate of evaluation point
%    y = y-coordinate of evaluation point
%  XY = nodal coordinates=[x y]
%OUTPUT;
%    s = s-coordinate of evaluation point
%    t = t-coordinate of evaluation point


function [s,t]=ISO_xy_to_st(x,y,XY)

NN=length(XY(:,1));


syms s t

N4 =[ 1/4-1/4*s-1/4*t+1/4*s*t
 1/4+1/4*s-1/4*t-1/4*s*t
 1/4+1/4*s+1/4*t+1/4*s*t
 1/4-1/4*s+1/4*t-1/4*s*t];

 N8 =[ -1/4+1/4*s*t+1/4*s^2+1/4*t^2-1/4*s^2*t-1/4*s*t^2
 -1/4-1/4*s*t+1/4*s^2+1/4*t^2-1/4*s^2*t+1/4*s*t^2
 -1/4+1/4*s*t+1/4*s^2+1/4*t^2+1/4*s^2*t+1/4*s*t^2
 -1/4-1/4*s*t+1/4*s^2+1/4*t^2+1/4*s^2*t-1/4*s*t^2
                      1/2-1/2*t-1/2*s^2+1/2*s^2*t
                      1/2+1/2*s-1/2*t^2-1/2*s*t^2
                      1/2+1/2*t-1/2*s^2-1/2*s^2*t
                      1/2-1/2*s-1/2*t^2+1/2*s*t^2];

                  % 

if NN==4;N8x=N4.';end
if NN==8;N8x=N8.';end

                  
                  fs=N8x*XY(:,1)
ft=N8x*XY(:,2);

EQ1=x-fs;
EQ2=y-ft;


SOL=(solve(EQ1,EQ2,s,t));
ss=double(SOL.s);
tt=double(SOL.t);
 I=find(ss>=-1 & ss<=1 & tt>=-1 & tt<=1);

s=ss(I);t=tt(I);


 
 

function  [XYData,EData]= M_Mesh_Q8_Function(XYData,Region)


%% Plot model regions
ITEXT=2;
figure;SeeRecModel(XYData,Region,ITEXT)
%disp('Line 7')
ITEXT=1;hold on;SeeRecModel(XYData,Region,ITEXT)
title('\bfRegions for mesh generation use Q8 elements')

X=XYData(:,2);Y=XYData(:,3);
%% Number of elements nd nodes per region
% NRegion=prod(Region(:,6:7)') % Number of nodes in region
% NN=cumsum(NRegion);
% ERegion=prod([(Region(:,6)-1) (Region(:,7)-1)]'); % Number of nodes in region
% NEN=cumsum(ERegion);

nn1=0;ne1=0;
                        %disp('Line 14 =================================')
TNodeXY=[]; % with duplicate coordinates
tALL=[];
%% Mesh all regions
for iR=1:length(Region(:,1));
Ns=Region(iR,6);S=linspace(-1,1,Ns);
Nt=Region(iR,7);T=linspace(-1,1,Nt);
%figure % Plot domain
IPLOT=0;
%[p,t,Elec]=RECT_MeshN(S,T,IPLOT);
[p,t,Elec,NodeXY]=RECT_8_MeshN(S,T,IPLOT);
%disp('Line 29')
%p=p
NXY=NodeXY;
p=NXY(:,2:3);


t=Elec(:,2:9);   %%%% 8-nodes for Q8 element

% size(t)
% size(Elec)
% t=t
% Elec=Elec
% pause
nn2=length(p(:,1))+nn1;
NN(iR)=nn2;
ne2=length(t(:,1))+ne1;
NEN(iR)=ne2;
ne1=ne2;nn1=nn2;
NRegion(iR)=length(p(:,1));
ERegion(iR)=length(t(:,1));
 %                       disp('Line 48 =================================') %,pause

%% Convert (s,t) into (x,y)
JABCD=Region(iR,2:5);
XE=X(JABCD);
YE=Y(JABCD);
XY=[XE YE];
ss=p(:,1);tt=p(:,2);
NodeXY=[];
for is=1:length(ss)
[x,y]=ISO_st_to_xy(ss(is),tt(is),XY);
NodeXY=[NodeXY;x y];
end
II=1:NRegion(iR);
if iR>1
II=(1:NRegion(iR))+NN(iR-1);
end
NodeTemp=[II' NodeXY];
TNodeXY=[TNodeXY;NodeTemp]; % with duplicate coordinates
if iR==1;nadd=0;end
if iR>1;nadd=NN(iR-1);end
tALL=[tALL;t+nadd];         % with node number to change
end
%disp('Line 71')
%% Define ElecTemp
IE=(1:length(tALL(:,1)))';
ElecTemp=[IE tALL];
disp('Line 585'),pause


%% Get rid of duplicated nodes
%% M_Mesh_Development_2.m
%% Get rid of duplicated nodes
XYT=TNodeXY;
format compact
XYKeep=[XYT(1:NN(1),[2 3])]; % keep all nodes in Region I
ILast=NN(1);           % Last node No.
NodeNew=(1:NN(1))';
for in=NN(1)+1:max(NN)
%for in=19:21
    XX=XYKeep(:,1);
    YY=XYKeep(:,2);
nK=length(XX);
   
     Xn=XYT(in,2);
     Yn=XYT(in,3);
  
     DX=(abs(Xn-XX));
     DY=(abs(Yn-YY));
    
     NI=find(DX+DY<1e-6);

     if length(NI)==0;XYKeep=[XYKeep;Xn Yn];
         ILast=ILast+1; NodeNew(in)=ILast;end
        
     if length(NI)>0;NodeNew(in)=NI;
     end
 end

%% Filtered nodes
IN=(1:length(XYKeep(:,1)))';
NodeXYNew=[IN XYKeep];
%% Change element connectivity data
% size(tALL)
J=tALL(:,1);JA=NodeNew(J);
  J=tALL(:,2);JB=NodeNew(J);
   J=tALL(:,3);JC=NodeNew(J);
    J=tALL(:,4);JD=NodeNew(J);
    J=tALL(:,5);
%disp('Line 117')
    NNN=NodeNew;
        JE=NodeNew(J);
    J=tALL(:,6);JF=NodeNew(J);
    J=tALL(:,7);JG=NodeNew(J);
    J=tALL(:,8);JH=NodeNew(J);
%    NEN,pause
    IE=(1:max(NEN))';
    ELECNew=[IE JA JB JC JD JE JF JG JH];
%% Plot mesh
XYData=NodeXYNew;
EData=ELECNew;
% %% write node number:
% ITEXT=1;
% SeeRecModel(XYData,EData,ITEXT)

%% write element number:
ITEXT=2; if max(NN)>30;ITEXT=0;end

figure
SeeRecModel(XYData,EData,ITEXT)
%% --------------------  end



%RECT_Mesh.m : Generate rectangular mesh
%function [p,t,Elec,X,Y]=RECT_Mesh(x1,x2,nx,y1,y2,ny,IPLOT)
function [p,t,Elec,X,Y]=RECT_MeshN(x,y,IPLOT)
% %Input
%     x1=0;x2=5;nx=4;
%     y1=0;y2=10;ny=6;

%     x=linspace(x1,x2,nx);
%     y=linspace(y1,y2,ny);
    nx=length(x);
    ny=length(y);
    
    
    
    [X,Y]=meshgrid(x,y);
XX=X;YY=Y;
%Element connectivity
        Elec=[];
    for IC=1:nx-1;

        I0  =(IC-1)*(ny-1);
        IE=I0+(1:ny-1)';

        JA0=(IC-1)*(ny);

        JA=(1:ny-1)'+JA0;
        JB=JA+ny;
        JC=JB+1;
        JD=JA+1;

Elec=[Elec;IE JA JB JC JD];
end
p=[X(:) Y(:)];
t=Elec(:,2:5);
                    tp=[t t(:,1)]; % repeat the first node, 7-19-12

NElem=length(Elec(:,1));
if IPLOT>0
    figure,plot(X(:),Y(:),'ro')
    text(X(:),Y(:),int2str((1:prod(size(X)))'))
hold on
end
for ie=1:NElem
    tN=tp(ie,:);
    xye=sum(p(tN,:));
    xy=sum(p(Elec(ie,2:5),:))/4;
    if IPLOT==1
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
        xye=p(tN,:);
    end
        
    if IPLOT>0
    plot(xye(:,1),xye(:,2),'b','linewidth',2)
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
end

end



% Modified from cc.m
function [p,t,Elec,X,Y]=Tri_Mesh_Cross2_Ring(x1,x2,NDx,y1,y2,NDy,IPLOT)
%% purpose: Mesh regions a ring type region bounded by 2 concentrated
%% circles and 2 radial lines
%% Special case : a disk defined by 2 concentric circles
% type  rect_mesh.m
% 
% %RECT_Mesh.m : Generate rectangular mesh
% function [p,t,Elec,X,Y]=RECT_Mesh(x1,x2,nx,y1,y2,ny,IPLOT)
% %Input
    yy=[5 6.5 7.5  9 10]; % user specified y 
%% Input
% NDx=8;NDy=2;
%     x1=0;x2=360; % angles
%     y1=5;y2=10;  % radius 
    ny=2*NDy+1;
nx=2*NDx+1;
    x=linspace(x1,x2,nx);
    y=linspace(y1,y2,ny);
    
        
    y=yy;
    
    [ANG,R]=meshgrid(x,y);
    
    X=R.*cosd(ANG);
    Y=R.*sind(ANG);
    
    
%% Dertermine center node number
ic=0;
t=[];
for ix=1:NDx
    for iy=1:NDy
        ic=ic+1;
        
        Na=2*(ix)*ny-ny;
        NC(ix,iy)=Na+(iy-1)*2+2;
        CN=Na+(iy-1)*2+2;
        CC(ic)=CN;
        A1=CN-(1+ny);
        A2=A1+1;
        A3=A2+1;
        
        A4=CN-1;
        A6=CN+1;
         A8=CN+ny;
         A7=A8-1;A9=A8+1;
         
         EL=[ CN A2 A1
             CN A1 A4
             CN A4 A7
             CN A7 A8
             CN A8 A9
             CN  A9 A6
             CN  A6 A3
             CN A3 A2];
     t=[t;EL];    
    end
end

        
        
        
        
    
    
    
    
    
NE=length(t(:,1));
IE=(1:NE)';
Elec=[IE t];

p=[X(:) Y(:)];
t=Elec(:,2:4);

NElem=length(Elec(:,1));
    figure,plot(X(:),Y(:),'ro')
    text(X(:),Y(:),int2str((1:prod(size(X)))'))
hold on

for ie=1:NElem
    tN=t(ie,:);
    xye=sum(p(tN,:));
    xy=sum(p(Elec(ie,2:4),:))/3;
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
        xye=p(tN,:);
 
        
%    if IPLOT>0
    plot(xye(:,1),xye(:,2),'b','linewidth',2)
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
  %end

end
ICC=CC;
plot(X(ICC),Y(ICC),'bs','linewidth',3)
% 
% [XYData,Elec,XX,YY]=Tri_Mesh_Cross2(x1,x2,NDx,y1,y2,NDy,IPLOT)
% title('\bfBy function')
if x2-x1==360
%%  M_ReplaceNode.m
NLast=nx*ny;
i1=NLast-ny+1;
NOLD=i1:NLast;
t0=t;
for in=1:ny
    for ic=1:3
        II=find(NOLD(in)==t(:,ic));
        t(II,ic)=II.^0*in;
    end
end
end


function  [XYData,EData]= M_Mesh_T3_Function(XYData,Region) % same as M_Mesh_Tri_function

%% Plot model regions
ITEXT=2;
figure;SeeRecModel(XYData,Region,ITEXT)
ITEXT=1;hold on;SeeRecModel(XYData,Region,ITEXT)
title('\bfRegions for triangular mesh generation')
X=XYData(:,2);Y=XYData(:,3);
%% Number of elements nd nodes per region
NRegion=prod(Region(:,6:7)'); % Number of nodes in region
NN=cumsum(NRegion);
ERegion=2*prod([(Region(:,6)-1) (Region(:,7)-1)]'); % Number of nodes in region
NEN=cumsum(ERegion);

TNodeXY=[]; % with duplicate coordinates
tALL=[];
%% Mesh region No.1

for iR=1:length(Region(:,1));
Ns=Region(iR,6);S=linspace(-1,1,Ns);
Nt=Region(iR,7);T=linspace(-1,1,Nt);
%figure % Plot domain
IPLOT=0;
% [p,t,Elec]=RECT_MeshN(S,T,IPLOT);
[p,t,Elec]=Tri_Mesh1_RECTN(S,T,IPLOT);
% disp('Line 20')
% pause
%% Convert (s,t) into (x,y)
JABCD=Region(iR,2:5);  %%%
XE=X(JABCD);
YE=Y(JABCD);
XY=[XE YE];
ss=p(:,1);tt=p(:,2);
NodeXY=[];
for is=1:length(ss)
[x,y]=ISO_st_to_xy(ss(is),tt(is),XY);
NodeXY=[NodeXY;x y];
end
II=1:NRegion(iR);
if iR>1
II=(1:NRegion(iR))+NN(iR-1);
end
NodeTemp=[II' NodeXY];
TNodeXY=[TNodeXY;NodeTemp]; % with duplicate coordinates
if iR==1;nadd=0;end
if iR>1;nadd=NN(iR-1);end
tALL=[tALL;t+nadd];         % with node number to change
end

%% Define ElecTemp
IE=(1:NEN(length(Region(:,1))))';
ElecTemp=[IE tALL];


%% Get rid of duplicated nodes
%% M_Mesh_Development_2.m
%% Get rid of duplicated nodes
XYT=TNodeXY;
format compact
XYKeep=[XYT(1:NN(1),[2 3])]; % keep all nodes in Region I
ILast=NN(1);           % Last node No.
NodeNew=(1:NN(1))';
for in=NN(1)+1:max(NN)
%for in=19:21
    XX=XYKeep(:,1);
    YY=XYKeep(:,2);
nK=length(XX);
   
     Xn=XYT(in,2);
     Yn=XYT(in,3);
  
     DX=(abs(Xn-XX));
     DY=(abs(Yn-YY));
    
     NI=find(DX+DY<1e-6);

     if length(NI)==0;XYKeep=[XYKeep;Xn Yn];
         ILast=ILast+1; NodeNew(in)=ILast;end
        
     if length(NI)>0;NodeNew(in)=NI;
     end
 end

%% Filtered nodes
IN=(1:length(XYKeep(:,1)))';
NodeXYNew=[IN XYKeep];
%% Change element connectivity data

J=tALL(:,1);JA=NodeNew(J);
  J=tALL(:,2);JB=NodeNew(J);
   J=tALL(:,3);JC=NodeNew(J);
%    J=tALL(:,4);JD=NodeNew(J);
    IE=(1:max(NEN))';
    ELECNew=[IE JA JB JC ];
%% Plot mesh
XYData=NodeXYNew;
EData=ELECNew;
% %% write node number:
% ITEXT=1;
% SeeRecModel(XYData,EData,ITEXT)
% 
% %% write element number:
 ITEXT=2;
 if max(NN)>30;ITEXT=0;end
% SeeRecModel(XYData,EData,ITEXT)
figure
SeeTriModel(XYData,EData,ITEXT);
%% --------------------  end



%diary M_EX13.m
%diary M_EX13.m For Example 5.1
%% To do:
%% 1. Generalize Logan9_2.m to an axymmetric program
%% 2. Use mesh generator to generate mesh for example 5.1
%% 3. Write a function M_EX13_General to solve the same problem 
%% with various meshes to show convergence of solutions

%% Plot Rectangulare model
function SeeRecModel(XYData,EData,ITEXT) %%%%%%%%%%%%%%%%%%
nin=nargin;
if nin==2;ITEXT=1;end

    I=XYData(:,1);
    NElem=length(EData(:,1));
X=XYData(I,2);
Y=XYData(I,3);
%   figure
plot(X(:),Y(:),'ro')
if ITEXT==1
     text(X(:),Y(:),int2str((1:prod(size(X)))'),'FontSize',20)
end
hold on
p=[X Y];
tp=EData(:,2:5);
 for ie=1:NElem
     tN=tp(ie,:); tN1=[tN tN(1)];%  , ITEXT,pause % repeat first node
     xye=(p(tN1,:));
     xy=sum(p(tN,:))/4;
     if ITEXT==2
     text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
         %xye=p(tN,:);
     end         
     plot(xye(:,1),xye(:,2),'b','linewidth',2)
% if ITEXT>0
%      text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',20)
% end
 end

% Same as : function M_Tri_Plot_Model(TIT,XYData,Elec,OPT)

%% In : C:\MATLAB701\work\2012-Summer-A\00-FE-Book-Work\00-FE-Book-BP-added\FirstDraft\Ch4-Aug-2012\Ch4-Codes\Tri-Codes

% Expand ranges
VV=axis;
RX=VV(2)-VV(1);RY=VV(4)-VV(3);
RR=max([RX RY]);
a=0.1;
V=VV+[-a a -a a]*RR; axis(V)
   

% %%  Tri_Mesh_Cross2_RingN.m
% type Tri_Mesh_Cross2_Ring
% 
% Modified from cc.m
function [p,t,Elec,X,Y]=Tri_Mesh_Cross2_RingN(x1,x2,NDx,y1,y2,NDy,IPLOT)
%% purpose: Mesh regions a ring type region bounded by 2 concentrated
%% circles and 2 radial lines
%% Special case : a disk defined by 2 concentric circles
% type  rect_mesh.m
% 
% %RECT_Mesh.m : Generate rectangular mesh
% function [p,t,Elec,X,Y]=RECT_Mesh(x1,x2,nx,y1,y2,ny,IPLOT)
% %Input
    yy=[5 6.5 7.5  9 10]; % user specified y 
%% Input
% NDx=8;NDy=2;
%     x1=0;x2=360; % angles
%     y1=5;y2=10;  % radius 
%     ny=2*NDy+1;
% nx=2*NDx+1;
%     x=linspace(x1,x2,nx);
%     y=linspace(y1,y2,ny);
    
        
    y=yy;
    
    [ANG,R]=meshgrid(x,y);
    
    X=R.*cosd(ANG);
    Y=R.*sind(ANG);
    
    
%% Dertermine center node number
ic=0;
t=[];
for ix=1:NDx
    for iy=1:NDy
        ic=ic+1;
        
        Na=2*(ix)*ny-ny;
        NC(ix,iy)=Na+(iy-1)*2+2;
        CN=Na+(iy-1)*2+2;
        CC(ic)=CN;
        A1=CN-(1+ny);
        A2=A1+1;
        A3=A2+1;
        
        A4=CN-1;
        A6=CN+1;
         A8=CN+ny;
         A7=A8-1;A9=A8+1;
         
         EL=[ CN A2 A1
             CN A1 A4
             CN A4 A7
             CN A7 A8
             CN A8 A9
             CN  A9 A6
             CN  A6 A3
             CN A3 A2];
     t=[t;EL];    
    end
end

        
        
        
        
    
    
    
    
    
NE=length(t(:,1));
IE=(1:NE)';
Elec=[IE t];

p=[X(:) Y(:)];
t=Elec(:,2:4);

NElem=length(Elec(:,1));
    figure,plot(X(:),Y(:),'ro')
    text(X(:),Y(:),int2str((1:prod(size(X)))'))
hold on

for ie=1:NElem
    tN=t(ie,:);
    xye=sum(p(tN,:));
    xy=sum(p(Elec(ie,2:4),:))/3;
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
        xye=p(tN,:);
 
        
%    if IPLOT>0
    plot(xye(:,1),xye(:,2),'b','linewidth',2)
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
  %end

end
ICC=CC;
plot(X(ICC),Y(ICC),'bs','linewidth',3)
% 
% [XYData,Elec,XX,YY]=Tri_Mesh_Cross2(x1,x2,NDx,y1,y2,NDy,IPLOT)
% title('\bfBy function')
if x2-x1==360
%%  M_ReplaceNode.m
NLast=nx*ny;
i1=NLast-ny+1;
NOLD=i1:NLast;
t0=t;
for in=1:ny
    for ic=1:3
        II=find(NOLD(in)==t(:,ic));
        t(II,ic)=II.^0*in;
    end
end
end




%% function  [N,Ns,Nt]=N_ISO48(s,t,NNode)
%% Purpose: 
%   Compute Shape function matrix for 4 and 8 node element
%% Usage: 
%       [N,Ns,Nt]=N_ISO48(s,t,NNode)
%% Input:  
%    s = s-coordinate of evaluation point
%    t = t-coordinate of evaluation point
% NNode= 4, for 4-noded element(default if not input)
%      =8, for 8-noded element
%OUTPUT;
%   N = shape function matrix
%   Ns= dN/ds (1X NNode) matrix
%   Nt= dN/dt (1X NNode) matrix
%% B,P, Wang
function  [N,Ns,Nt]=N_ISO48(s,t,NNode)
%Modified from  function [N,Ns,Nt]=K48_st(s,t,NNode)

if nargin==2;NNode=4; end % default= 4 node element

if NNode==4
N4 =[ 1/4-1/4*s-1/4*t+1/4*s*t
 1/4+1/4*s-1/4*t-1/4*s*t
 1/4+1/4*s+1/4*t+1/4*s*t
 1/4-1/4*s+1/4*t-1/4*s*t];
 P=[1-s^2  1-t^2];

 N=N4.';
                        %N4s=diff(N4,s)
 N4s =[ -1/4+1/4*t
  1/4-1/4*t
  1/4+1/4*t
 -1/4-1/4*t];
 %N4t=diff(N4,t)
 Ns=N4s.';
 N4t =[ -1/4+1/4*s
 -1/4-1/4*s
  1/4+1/4*s
  1/4-1/4*s];
 Nt=N4t.';
end
if NNode==8

 %N8t=diff(N8,t)
 N8t =[   1/4*s+1/2*t-1/4*s^2-1/2*s*t
 -1/4*s+1/2*t-1/4*s^2+1/2*s*t
  1/4*s+1/2*t+1/4*s^2+1/2*s*t
 -1/4*s+1/2*t+1/4*s^2-1/2*s*t
                 -1/2+1/2*s^2
                       -t-s*t
                  1/2-1/2*s^2
                       -t+s*t];
 Nt=N8t.';
 % N8s=diff(N8,s)
 N8s =[  1/4*t+1/2*s-1/2*s*t-1/4*t^2
 -1/4*t+1/2*s-1/2*s*t+1/4*t^2
  1/4*t+1/2*s+1/2*s*t+1/4*t^2
 -1/4*t+1/2*s+1/2*s*t-1/4*t^2
                       -s+s*t
                  1/2-1/2*t^2
                       -s-s*t
                 -1/2+1/2*t^2];
 Ns=N8s.';
 N8 =[ -1/4+1/4*s*t+1/4*s^2+1/4*t^2-1/4*s^2*t-1/4*s*t^2
 -1/4-1/4*s*t+1/4*s^2+1/4*t^2-1/4*s^2*t+1/4*s*t^2
 -1/4+1/4*s*t+1/4*s^2+1/4*t^2+1/4*s^2*t+1/4*s*t^2
 -1/4-1/4*s*t+1/4*s^2+1/4*t^2+1/4*s^2*t-1/4*s*t^2
                      1/2-1/2*t-1/2*s^2+1/2*s^2*t
                      1/2+1/2*s-1/2*t^2-1/2*s*t^2
                      1/2+1/2*t-1/2*s^2-1/2*s^2*t
                      1/2-1/2*s-1/2*t^2+1/2*s*t^2];
 N=N8.';
end
%---------------- end of N_ISO48.m




%diary M_EX13.m
%diary M_EX13.m For Example 5.1
%% To do:
%% 1. Generalize Logan9_2.m to an axymmetric program
%% 2. Use mesh generator to generate mesh for example 5.1
%% 3. Write a function M_EX13_General to solve the same problem 
%% with various meshes to show convergence of solutions

%% Plot triangle model
function SeeTriModel(XYData,EData,ITEXT) %%%%%%%%%%%%%%%%%%
nin=nargin;
if nin==2;ITEXT=1;end

    I=XYData(:,1);
    NElem=length(EData(:,1));
X=XYData(I,2);
Y=XYData(I,3);
%   figure
plot(X(:),Y(:),'ro')
if ITEXT>0
     text(X(:),Y(:),int2str((1:prod(size(X)))'),'FontSize',20)
end
hold on
p=[X Y];
tp=EData(:,2:4);
 for ie=1:NElem
     tN=tp(ie,:); tN1=[tN tN(1)];%  , ITEXT,pause % repeat first node
     xye=(p(tN1,:));
     xy=sum(p(tN,:))/3;
%%     text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
         %xye=p(tN,:);
         
     plot(xye(:,1),xye(:,2),'b','linewidth',2)
if ITEXT>0
     text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',20)
end
 end

% Same as : function M_Tri_Plot_Model(TIT,XYData,Elec,OPT)

%% In : C:\MATLAB701\work\2012-Summer-A\00-FE-Book-Work\00-FE-Book-BP-added\FirstDraft\Ch4-Aug-2012\Ch4-Codes\Tri-Codes

% Expand ranges
VV=axis;
RX=VV(2)-VV(1);RY=VV(4)-VV(3);
RR=max([RX RY]);
a=0.1;
V=VV+[-a a -a a]*RR; axis(V)
   

%diary RECT_4_MeshN.m



%% function [p,t,Elec,NodeXY]=RECT_4_Mesh(x1,x2,nx,y1,y2,ny,IPLOT);
%% Purpose:
%   Generate uniform rectangular mesh for 4-noded rectangular elements in a rectangular domain
%% Input
% x1,x2,nx   = minimum x, maximun x, number of points
% y1, y2, ny = minimum y, maximun y, number of points
%   IPLOT =1; % plot medh
%         =0, do not plot
%% Output
%   p  = [x y] nodal coordinates
%   t  = [JA JB JC JD], element connectivity
%  Elec= [ElenentNo JA JB JC JD ], element connectivity
%  NodeXY= [NodeNo X  Y], element connectivity
%
%% B.P. Wang, July 19,2012
% In : C:\MATLAB701\work\2012-Summer-A\00-FE-Book-Work\00-FE-Book-BP-added\Ch7-Mesh

%function [p,t,Elec,NodeXY]=RECT_4_Mesh(x1,x2,nx,y1,y2,ny,IPLOT);
function [p,t,Elec,NodeXY]=RECT_4_MeshN(x,y,IPLOT);

% Developed in REC8_Mesh_Development.m
[p,t,Elec]=RECT_MeshN(x,y,IPLOT);
% aa=0.2;AxisExpand(aa)
% %% ===================================================
% %% Add vertical nodes
% DY=(y2-y1)/(ny-1);
%     x=linspace(x1,x2,nx);
%     y=linspace(y1+DY/2,y2-DY/2,ny-1);
%     [X,Y]=meshgrid(x,y);
% px=X(:);py=Y(:);
% pp=[p;px py]
% plot(px,py,'ro','markersize',12)
% nv=nx*(ny-1);
% np=nx*ny+(1:nv);
% text(px,py,int2str(np'),'color','r','Fontsize',15)
% %% Add horizonztal nodes
% DX=(x2-x1)/(nx-1);
%     y=linspace(y1,y2,ny);
%     x=linspace(x1+DX/2,x2-DX/2,nx-1);
%     [X,Y]=meshgrid(x,y);
% px=X(:);py=Y(:);
% pp=[pp;px py]
% plot(px,py,'rs','markersize',12)
% nh=ny*(nx-1);
% np=nx*ny+nv+(1:nh);
% text(px,py,int2str(np'),'color','r','Fontsize',15)
% 
% %% Nodes H and F
% nn=nx*ny;
% nnx=nx*(ny-1);
% nny=ny*(nx-1);
% HH=nn+(1:(nx-1)*(ny-1))';
% FF=HH+(ny-1);
% [HH FF]
% 
% %% Nodes E
% EE=[];
% for ix=1:nx-1
%     e1=(nn+nnx)+(ix-1)*ny+(1:(ny-1))';
%     EE=[EE;e1];
% end
% %% Nodes G
% GG=EE+1;
% t58=[EE FF GG HH];
% t8=[t t58];
 NE=length(t(:,1));
%% Output
Elec=[(1:NE)' t ones(NE,2)];
II=(1:length(p(:,1)))';
NodeXY=[II p];




%% function  [p,t,Elec,X,Y]=Tri_Mesh1_RECTN(x,y,IPLOT)
%% Purpose:
%      Generate triangular mesh rectangular domain ('+45')
%% Usage:
%      [p,t,Elec,X,Y]=Tri_Mesh1_RECTN(x,y,IPLOT)
%%Input
%     x = x-xoordinates for mesh points
%     y = y-xoordinates for mesh points

function  [p,t,Elec,X,Y]=Tri_Mesh1_RECTN(x,y,IPLOT)

% %RECT_Mesh.m : Generate rectangular mesh
% function [p,t,Elec]=RECT_Mesh(x1,x2,nx,y1,y2,ny,IPLOT)
% IPLOT=0, do not plot
%      =1, plot model
%     x=linspace(x1,x2,nx);
%     y=linspace(y1,y2,ny);
    [X,Y]=meshgrid(x,y);
nx=length(x);
ny=length(y);

%Element connectivity
        Elec=[];
    for IC=1:nx-1;

        I0  =2*(IC-1)*(ny-1);
        IE=I0+(1:2*(ny-1))';

        JA0=(IC-1)*(ny);

        JA=(1:ny-1)'+JA0;
        JB=JA+ny;
        JC=JB+1;
        JD=JA+1;
  TT=[JA JB JC
      JA JC JD];
  Elec=[Elec;IE TT ];
end
p=[X(:) Y(:)];
t=Elec(:,2:4);

NElem=length(Elec(:,1));
%IPLOT=1
if IPLOT>0
figure,plot(X(:),Y(:),'ro')
    text(X(:),Y(:),int2str((1:prod(size(X)))'))
hold on
end
for ie=1:NElem
    tN=[t(ie,:) t(ie,1)];
    xye=p(tN,:);
    xy=sum(p(Elec(ie,2:4),:))/3;
    if IPLOT>0
    plot(xye(:,1),xye(:,2),'b','linewidth',2)
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
end
end




%diary RECT_8_MeshN.m
%type RECT_8_Mesh

%% function [p,t,Elec,NodeXY]=RECT_8_MeshN(x1,x2,nx,y1,y2,ny,IPLOT);
%% Purpose:
%   Generate uniform rectangular mesh for 8-noded rectangular elements in a rectangular domain
%% Input
% x1,x2,nx   = minimum x, maximun x, number of points
% y1, y2, ny = minimum y, maximun y, number of points
%   IPLOT =1; % plot medh
%         =0, do not plot
%% Output
%   p  = [x y] nodal coordinates
%   t  = [JA JB JC JD], element connectivity
%  Elec= [ElenentNo JA JB JC JD JE JF JG JH], element connectivity
%  NodeXY= [NodeNo X  Y], element connectivity
%
%% B.P. Wang, July 19,2012

function [p,t,Elec,NodeXY]=RECT_8_MeshN(x,y,IPLOT);

% Developed in REC8_Mesh_Development.m
[p,t,Elec]=RECT_MeshN(x,y,IPLOT);
%disp('Line 24'),pause

% % Expand ranges
% VV=axis;
% RX=VV(2)-VV(1);RY=VV(4)-VV(3);
% RR=max([RX RY]);
% a=0.1;
% V=VV+[-a a -a a]*RR; axis(V)
% 

%aa=0.2;AxisExpand(aa)
%% ===================================================
%% Add vertical nodes
%% Unidorm
% DY=(y2-y1)/(ny-1);  %% to be done
%     x=linspace(x1,x2,nx);
%     y=linspace(y1+DY/2,y2-DY/2,ny-1);
%     [X,Y]=meshgrid(x,y);
% px=X(:);py=Y(:);
% pp=[p;px py]
%% Non-Unidorm
ny=length(y);
nx=length(x);
DY=y(2:ny)-y(1:(ny-1));  %% to be done

ya=y(1:(ny-1))+DY/2;
     [X,Y]=meshgrid(x,ya);
 px=X(:);py=Y(:);
 pp=[p;px py];



if IPLOT>0
    plot(px,py,'ro','markersize',12)
nv=nx*(ny-1);
np=nx*ny+(1:nv);
text(px,py,int2str(np'),'color','r','Fontsize',15)
end
%% Add horizonztal nodes
%% Uniform
% DX=(x2-x1)/(nx-1);
%     y=linspace(y1,y2,ny);
%     x=linspace(x1+DX/2,x2-DX/2,nx-1);
%     [X,Y]=meshgrid(x,y);
%% Non-Uniform

 DX=x(2:nx)-x(1:(nx-1));

 x=x(1:nx-1)+DX/2;
 [X,Y]=meshgrid(x,y);
px=X(:);py=Y(:);
pp=[pp;px py];
if IPLOT>0
    plot(px,py,'rs','markersize',12)
nh=ny*(nx-1);
np=nx*ny+nv+(1:nh);
text(px,py,int2str(np'),'color','r','Fontsize',15)
end
%% Nodes H and F
nn=nx*ny;
nnx=nx*(ny-1);
nny=ny*(nx-1);
HH=nn+(1:(nx-1)*(ny-1))';
FF=HH+(ny-1);
%[HH FF]

%% Nodes E
EE=[];
for ix=1:nx-1
    e1=(nn+nnx)+(ix-1)*ny+(1:(ny-1))';
    EE=[EE;e1];
end
%% Nodes G
GG=EE+1;
t58=[EE FF GG HH];
t8=[t t58];
NE=length(t8(:,1));
%% Output
Elec=[(1:NE)' t8 ones(NE,2)];
II=(1:length(pp(:,1)))';
NodeXY=[II pp];




% %diary Tri_Mesh2_RECTN.m
% type  Tri_Mesh2_RECT
% 

function  [p,t,Elec,X,Y]=Tri_Mesh2_RECTN(x,y,IPLOT)
% Generate triangular mesh rectangular domain ('-45')
%               : a special case
% %Input
%     x1=0;x2=5;nx=4;
%     y1=0;y2=10;ny=6;

% %RECT_Mesh.m : Generate rectangular mesh
% function [p,t,Elec]=RECT_Mesh(x1,x2,nx,y1,y2,ny,IPLOT)
% IPLOT=0, do not plot
% %      =1, plot model
%     x=linspace(x1,x2,nx);
%     y=linspace(y1,y2,ny);
nx=length(x);ny=length(y);
    [X,Y]=meshgrid(x,y);

%Element connectivity
        Elec=[];
    for IC=1:nx-1;

        I0  =2*(IC-1)*(ny-1);
        IE=I0+(1:2*(ny-1))';

        JA0=(IC-1)*(ny);

        JA=(1:ny-1)'+JA0;
        JB=JA+ny;
        JC=JB+1;
        JD=JA+1;
  TT=[JA JB JD
      JB JC JD];
  Elec=[Elec;IE TT ];
end
p=[X(:) Y(:)];
t=Elec(:,2:4);

NElem=length(Elec(:,1));
IPLOT=1
if IPLOT>0
figure,plot(X(:),Y(:),'ro')
    text(X(:),Y(:),int2str((1:prod(size(X)))'))
hold on
end
for ie=1:NElem
    tN=[t(ie,:) t(ie,1)];
    xye=p(tN,:);
    xy=sum(p(Elec(ie,2:4),:))/3;
% p(Elec(ie,2:4),:)
% ie=ie
  
    if IPLOT>0
    plot(xye(:,1),xye(:,2),'b','linewidth',2)
    text(xy(1),xy(2),int2str(ie),'Color','b','FontSize',12)
  %pause
end
end

% -----------------------
% See stand-alone version
% -----------------------
% function [NLine]=Find_Line_Node(KeyPoint ,XYData,XYData1) % May-12-2015
% %% Purpose: 
% %  Find nodes between key-points
% %% Input
% 
% % XYData =Nodes to define regions
% %       =[ Node  X    Y]
% 
% %% KeyPoint=[KP1 KP2]=ID of key points (nodes in CYData);
% %% XYData1 = data generated by M_Mesh2_A or B
% 
% XYK1=XYData(KeyPoint(1),[2 3]);
% XYK2=XYData(KeyPoint(2),[2 3]);
% XYK=[XYK1-XYK2];
%  LXY=norm(XYK);
% %% check by finding the distance by cross product
% IN=XYData1(:,1);
% 
% NNode=length(IN);
% ILine=[];  % nodes on this line
% for i=1:NNode
%     In=IN(i);
%     XY=XYData1(In,[2 3]);
%    
%     V1=[XYK1-XY 0];
%     V2=[XYK2-XY 0];
%     V12=cross(V1,V2);
% %     V1X=abs(sum(V1.*[XYK 0]))/LXY;
% %     V2X=abs(sum(V2.*[XYK 0]))/LXY;
%     V1X=abs(sum(V1.*[XY 0]))/LXY;  % 5-12-15
%     V2X=abs(sum(V2.*[XY 0]))/LXY;
%     if abs(V12(3))<1e-6 & (V1X+V2X)<=LXY; ILine=[ILine In]; end
% end
% disp(['The generated nodes between keypoints[ ',int2str(KeyPoint),' ] are:'])
%     disp(ILine)
% NLine=ILine;

%% end of all   sub-functions ================================ 

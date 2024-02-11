function ISO_drawArc(XY4)
%% ISO_drawArc.m
% The curve is a 270 degree arc of a circle
% in : C:\2014-FEM-AE-ME-5310\Daily-Notes\Nov-2014\ISO-TOOls-Nov-2014\XY-to-st-Feb17-2015
%% An example of computing (s,t) from (x,y)
%% BPW, Feb 18,2015

%% Draw a arc in x-y plane
% set circle origin and radius
%XC=0;YC=10;R=6;th=linspace(pi/2,2*pi,200);
[XC,YC]=ISO_st_to_xy(0.0,0.0,XY4);

% Est a radius inside quad as 0.25 min(a,b)
a=abs( ( max( XY4(:,1) ) -min( XY4(:,1) ) ) );
b=abs( ( max( XY4(:,2) ) -min( XY4(:,2) ) ) );
R=0.25*(min(a,b));

% Define a 270 degree arc
th=linspace(pi/2,2*pi,200);
xc=XC+R*cos(th);
yc=YC+R*sin(th);
subplot(1,2,1)
hold on
plot(xc,yc,'r','linewidth',2)

text(xc(1),yc(1),'\bfA')
text(xc(end),yc(end),'\bfB')
text(XC,YC-.82,'\bfC')

%% Draw arc in s-t plane
for i=1:length(th)
    Px=xc(i);Py=yc(i);
    % [Ps,Pt]=ISO_Points_xy_to_st_NoPlot(XY4,Px,Py,PLab);
    [Ps,Pt]=ISO_XY_to_st_Num(Px,Py,XY4);
    PS(i)=Ps;PT(i)=Pt;
end
subplot(1,2,2)
hold on
plot(PS,PT,'r','linewidth',2)
[PsC,PtC]=ISO_XY_to_st_Num(XC,YC,XY4);

text(PS(1),PT(1),'\bfA')
text(PS(end),PT(end),'\bfB')
text(PsC,PtC-.12,'\bfC')


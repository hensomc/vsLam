%function ISO_PlotQ4_Sub(XY) %  Plot Q4 element with (s,t) axis
%  Plot Q4 element in physical and natural coordinates
%Input:
%   XY= [x y] coordinates of a Q4 element
%--------------------------------------------------
%Modified from ISO_Stydy_1.m
%% B.P. Wang, Oct 20,2014

function ISO_PlotQ4(XY) %  Plot Q4 element in (x,y) & (s,t) systems

% draw quad in x-y and corner node IDs
XY4=XY;
xy=[XY4;XY4(1,:)];
subplot(1,2,1)
plot(xy(:,1),xy(:,2),'b','linewidth',3)
text(XY4(:,1)+.5,XY4(:,2),int2str((1:4)'),'color','r','fontsize',14)
axis equal;

%% Plot x-axis (transformed s-axis)
ss=[-1.2;1.2];tt=[0;0];
for i=1:2
    s=ss(i);t=tt(i);
    [x,y]=ISO_st_to_xy(s,t,XY4);
    XA(i,1)=x;
    YA(i,1)=y;
end
hold on;
subplot(1,2,1)
plot(XA,YA,'b--','linewidth',2)
text(XA(2,1),YA(2,1)+0.25,'\bfx','FontSize',14)

%% Plot y-axis (transformed t-axis)
tt=[-1.2;1.2];ss=[0;0];
for i=1:2
    s=ss(i);t=tt(i);
    [x,y]=ISO_st_to_xy(s,t,XY4);
    XA(i,1)=x;
    YA(i,1)=y;
end
hold on;
subplot(1,2,1)
plot(XA,YA,'b--','linewidth',2)
text(XA(2,1),YA(2,1)+0.25,'\bfy','FontSize',14)

%% Plot origin in x-y (transformed from s-t)
subplot(1,2,2)
s=0;t=0;NN=4;
N=N_ISO48(s,t,NN);
X1=s;
Y1=t;
hold on;
plot(X1,Y1,'ro','linewidth',3)
   
%% draw unit Quad in s-t and corner node IDs
XY4=[-1 -1; 1 -1; 1 1; -1 1];
xy=[XY4;XY4(1,:)];
plot(xy(:,1),xy(:,2),'r','linewidth',2)
text(XY4(:,1)+.15,XY4(:,2),int2str((1:4)'),'color','r','fontsize',14)
axis equal;

%% plot s-axis
s=-1.2;t=0;NN=4;
plot([s -s],[0 0],'r--','linewidth',1.5)
%text(-s,.1,'\bfs','FontSize',14)
text(-s,.1,'\bf\xi','FontSize',14)

%% plot t-axis
t=-1.2;s=0;NN=4;
plot([0 0],[t -t],'r--','linewidth',1.5)
%text(0,-t,'\bft','FontSize',14)
text(0,-t,'\bf\eta','FontSize',14)


 
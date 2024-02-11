%ISO_xy_to_st.m Given (x,y) compute (s,t)
%% function [s,t]=ISO_xy_to_st(x,y,XY)
%% Purpose:
%   Given (x,y) compute (s,t) in Q4 or Q8 element
%% Usage : 
%  [s,t]=ISO_xy_to_st(x,y,XY)
%% Input:  
%    x = x-coordinate of evaluation point
%    y = y-coordinate of evaluation point
%  XY = nodal coordinates=[x y]
%% OUTPUT;
%    s = s-coordinate of evaluation point
%    t = t-coordinate of evaluation point
%% Remarks: matlab numerical solution function fsolve  is used in this program.

%% B.P. Wang, Feb 18,2015
%% Test1.m
function [s,t]=ISO_XY_to_st_Num(x,y,XY)
% XY4=[             -20             0
%          40.00          8.00
%           20.00          25.00
%           -35.00          14.00];
% XC=0;YC=10;
ST0=[0;0];
% Compute start point as panel centroid
%ST0=[mean(XY(:,1)); mean(XY(:,2))]

options = optimset('Display','off');
st=fsolve(@ISO_xy_to_st_Fun,ST0,[options],x,y,XY);
s=st(1);t=st(2);

function [f]=ISO_xy_to_st_Fun(ST,x,y,XY)
% comment: make sure point on element edge is insode the element by a smalll amount

s=ST(1);t=ST(2);
N4 =[ 1/4-1/4*s-1/4*t+1/4*s*t
 1/4+1/4*s-1/4*t-1/4*s*t
 1/4+1/4*s+1/4*t+1/4*s*t
 1/4-1/4*s+1/4*t-1/4*s*t];

fs=N4'*XY(:,1);
ft=N4'*XY(:,2);

EQ1=x-fs;
EQ2=y-ft;
f=[EQ1;EQ2];
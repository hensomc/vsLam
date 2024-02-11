%% function [x,y]=ISO_st_to_xy(s,t,XY)
%% Purpose:
%   Given (s,t) compute (x,y) in Q4 or Q8 element
%% Usage:
%  [x,y]=ISO_st_to_xy(s,t,XY)
%% Input:  
%    s = s-coordinate of evaluation point
%    t = t-coordinate of evaluation point
%  XY = nodal coordinates=[x y]
%      XY = 4 X 2, for Q4 element
%         = 8 X 2  for Q8 element
%% OUTPUT;
%    x = x-coordinate of evaluation point
%    y = y-coordinate of evaluation point

%% B.P. Wang, Oct 12,2012
function [x,y]=ISO_st_to_xy(s,t,XY)

NN=length(XY(:,1));
    N=N_ISO48(s,t,NN);
    x=N*XY(:,1);
    y=N*XY(:,2);
   
   
%% function  [N,Ns,Nt]=N_ISO48(s,t,NNode)
%% Purpose:
%   Compute shape function (N)  matrix for 4 and 8 node quadrilateral elements 
%% Input:  
%    s = s-coordinate of evaluation point
%    t = t-coordinate of evaluation point
% NNode= 4, for 4-noded element(default if not input)
%      =8, for 8-noded element
%% OUTPUT;
%   N = shape function matrix
%   Ns= dN/ds (1X NNode) matrix
%   Nt= dN/dt (1X NNode) matrix

%% B.P. Wang, Oct 12,2012

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
%---------------- nd of N_ISO48.m



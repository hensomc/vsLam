function [F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_power(Nmax,xn,yn)
%% Purpose:  Generate F and its derivative, numerically, vectorized operation
% F=[1 x y x^2 x*y y^2 ...] = terms in the Pascal triangle that includes 
% combination of terms (1,x,x^2 ...x^Nmax)
%                and   (1,y,y^2 ...y^Nmax)
%% Input:
%   Nmax = maximum order of polynomial terms, in x and in y
%   xn = x-coordinates of evaluation points
%   yn = y-coordinates of evaluation points
%% Output:
%   F = [xn.^0  xn yn xn.^2 xn.*yn yn.^2 ....], a matrix (Nx X Nterm)
%        where Nx = number of evaluation points(length of vector xn)
%              Nterm=Nmax*(Nmax+1)/2, number of terms in polynomial 
%  dFx = dF/dx at evaluation points,a matrix (Nx X Nterm)
%  dFy = dF/dy at evaluation points,a matrix (Nx X Nterm)
%  dFxy = d^2F/dxdy at evaluation points,a matrix (Nx X Nterm)
%  dFxx = d^2F/dx^2 at evaluation points,a matrix (Nx X Nterm)
%  dFyy = d^2F/dy^2 at evaluation points,a matrix (Nx X Nterm)
%
%%  B.P. Wang, MAE UTA,May 7,2014
for i=1:length(xn)
FXY=[1];
    for N=1:Nmax
    x=xn(i);y=yn(i);
    J=0:N;I=N-J;
xx=x.^I;yy=y.^J;FXY=[FXY xx.*yy];end
FF(i,:)=FXY;
end
F=FF;

%For dF/dx
for i=1:length(xn)
    FXY=[0];
    x=xn(i);y=yn(i);
    for N=1:Nmax
        J=0:N;
        I=N-J;
        xx=I.*x.^(I-1);  % undefined when x=0 and I = 0
        yt=y;
        yy=yt.^J;
        FXY=[FXY xx.*yy];
    end
    dFx(i,:)=FXY;
end

%For dF/dy
for i=1:length(xn)
    FXY=[0];
    x=xn(i);y=yn(i);
    for N=1:Nmax
        J=0:N;
        I=N-J;
        xx=x.^I;
        yy=J.*y.^(J-1);  % undefined when y=0 and J= 0 ??
        FXY=[FXY xx.*yy];
    end
    dFy(i,:)=FXY;
end

%For dFxx
for i=1:length(xn)
    FXY=[0];
    x=xn(i);y=yn(i);
    for N=1:Nmax
        J=0:N;
        I=N-J;
        xx=I.*(I-1).*x.^(I-2);
        yy=y.^J;
        
        FXY=[FXY xx.*yy];end
    dFxx(i,:)=FXY;
end

%For dFxy
for i=1:length(xn)
    FXY=[0];
    x=xn(i);y=yn(i);
    for N=1:Nmax
        J=0:N;
        I=N-J;
         xx=I.*x.^(I-1);
        yy=J.*y.^(J-1);
        FXY=[FXY xx.*yy];end
    dFxy(i,:)=FXY;
end

%For dFyy
for i=1:length(xn)
    FXY=[0];
    x=xn(i);y=yn(i);
    for N=1:Nmax
        J=0:N;
        I=N-J;
         xx=x.^I;
        yy=J.*(J-1).*y.^(J-2);
        FXY=[FXY xx.*yy];end
    dFyy(i,:)=FXY;
end
%% ==============================  all done

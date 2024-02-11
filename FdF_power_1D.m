function [F,dFx,dFxx]=FdF_power_1D(Nmax,xn)
%% Purpose:  Generate F and its derivative, numerically, vectorized operation
% F=[1 x x^2 x^3 ...] = terms in the Pascal triangle that includes 
% combination of terms (1,x,x^2 ...x^Nmax)

%% Input:
%   Nmax = maximum order of polynomial terms, in x
%   xn = x-coordinates of evaluation points

%% Output:
%   F = [xn.^0  xn xn.^2 ....], a matrix (Nx X Nterm)
%        where Nx = number of evaluation points(length of vector xn)
%              Nterm=Nmax*(Nmax+1)/2, number of terms in polynomial 
%  dFx = dF/dx at evaluation points,a matrix (Nx X Nterm)
%  dFxx = d^2F/dx^2 at evaluation points,a matrix (Nx X Nterm)
%
%%  B.P. Wang, MAE UTA,May 7,2014
 
%For F
%disp('FdF_power_1D called');
for i=1:length(xn)
    FXY=[1];
%     for N=1:Nmax
%         x=xn(i);
%         J=0:N;
%         I=N-J;
%         xx=x.^I;FXY=[FXY xx];
%     end

    for I=1:Nmax
        x=xn(i);
        xx=x.^I;
        FXY=[FXY xx];
    end
    FF(i,:)=FXY;
end
F=FF;

%For dF/dx
for i=1:length(xn)
    FXY=[0 1];
    x=xn(i);
%     for N=1:Nmax
%         J=0:N;
%         I=N-J;
%         xx=I.*x.^(I-1);  % undefined when x=0 and I = 0
%         FXY=[FXY xx];
%     end
    for I=2:Nmax
        xx=I.*x.^(I-1);  % undefined when x=0 and I = 0
        FXY=[FXY xx];
    end
    dFx(i,:)=FXY;
end


%For dFxx
for i=1:length(xn)
    FXY=[0 0 2];
    x=xn(i);
%     for N=1:Nmax
%         J=0:N;
%         I=N-J;
%         xx=I.*(I-1).*x.^(I-2);        
%         FXY=[FXY xx];end
    for I=3:Nmax
        xx=I.*(I-1).*x.^(I-2);        
        FXY=[FXY xx];
    end
    dFxx(i,:)=FXY;
end

%% ==============================  all done

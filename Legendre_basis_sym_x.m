function [LEGFUN,LEGNx]=Legendre_basis_sym_x(NN)
% Generate NN-th order Legendre polynomials in 's' and 'x' domain
%% Input
%   NN = order
%% Output:
% LEGFUN= Legendre polynomial in s domain   % -1<=s<=1
% LEGNx=Legendre polynomial in x domain % 0<=x<=L
%% B.P. Wang, MAE UTA, Jan 27, 2015
% in : C:\2015\6311\Ritz-Paper
%% in C:\2015\6311\Ritz-Paper

% B.P. wang, MAE UTA, 2015
LEGFUN=[1];
syms x s L
for i=1:NN
    LEG=Legendre_Function(i);
    LEGFUN=[LEGFUN LEG];
end
%% Convert into function of x, x=0 to L
 LEGFUN= subs(LEGFUN,x,s);% Legendre polynomial in s domain   % -1<=s<=1
 LEGNx=subs(LEGFUN,s, (2*x/L-1) );%Legendre polynomial in x domain % 0<=x<=L

function [LEGN]=Legendre_Function(n)

%C:\2015\6311\Ritz-Beam-More-jan19-2015

syms x 
%for nn=1:n
Pn=0;
for i=1:(n+1);
    k=i-1;a=(-1)^k;b=Cnk_Fun(n,k);b=b^2;f1=(1+x)/2;f2=(1-x)/2;
    Pi=a*b*f1^(n-k)*f2^k;
    PLEG(i)=Pi;
    Pn=Pn+Pi;
end

LEGN=simplify(Pn);
function Cnk=Cnk_Fun(N,k)
FN=factorial(N);
FK=factorial(k);
FNMK=factorial(N-k);
Cnk=FN/(FK*FNMK);

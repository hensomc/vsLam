function [F,Fx,Fy,Fxy,Fxx,Fyy] = FdF_legendre3(nDeg,x,y)

F   = F_legendre_3(x,y);
Fx  = Fx_legendre_3(x,y);
Fy  = Fy_legendre_3(x,y);
Fxy = Fxy_legendre_3(x,y);
Fxx = Fxx_legendre_3(x,y);
Fyy = Fyy_legendre_3(x,y);


function F = F_legendre_3(x,y)
%F_LEGENDRE_3
%    F = F_LEGENDRE_3(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    14-Apr-2015 22:28:40

t2 = y.^2;
t3 = t2.*(3.0./2.0);
t4 = t3-1.0./2.0;
t5 = t2.*5.0;
t6 = t5-3.0;
t7 = x.^2;
t8 = t7.*(3.0./2.0);
t9 = t8-1.0./2.0;
t10 = t7.*5.0;
t11 = t10-3.0;
F = [1.0,y,t4,t6.*y.*(1.0./2.0),x,x.*y,t4.*x,t6.*x.*y.*(1.0./2.0),t9,t9.*y,t4.*t9,t6.*t9.*y.*(1.0./2.0),t11.*x.*(1.0./2.0),t11.*x.*y.*(1.0./2.0),t4.*t11.*x.*(1.0./2.0),t6.*t11.*x.*y.*(1.0./4.0)];


function Fx = Fx_legendre_3(x,y)
%FX_LEGENDRE_3
%    FX = FX_LEGENDRE_3(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    14-Apr-2015 22:28:40

t2 = y.^2;
t3 = t2.*(3.0./2.0);
t4 = t3-1.0./2.0;
t5 = t2.*5.0;
t6 = t5-3.0;
t7 = x.^2;
t8 = t7.*(1.5e1./2.0);
t9 = t8-3.0./2.0;
Fx = [0.0,0.0,0.0,0.0,1.0,y,t4,t6.*y.*(1.0./2.0),x.*3.0,x.*y.*3.0,t4.*x.*3.0,t6.*x.*y.*(3.0./2.0),t9,t9.*y,t4.*t9,t6.*t9.*y.*(1.0./2.0)];


function Fy = Fy_legendre_3(x,y)
%FY_LEGENDRE_3
%    FY = FY_LEGENDRE_3(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    14-Apr-2015 22:28:40

t2 = y.^2;
t3 = t2.*(1.5e1./2.0);
t4 = t3-3.0./2.0;
t5 = x.^2;
t6 = t5.*(3.0./2.0);
t7 = t6-1.0./2.0;
t8 = t5.*5.0;
t9 = t8-3.0;
Fy = [0.0,1.0,y.*3.0,t4,0.0,x,x.*y.*3.0,t4.*x,0.0,t7,t7.*y.*3.0,t4.*t7,0.0,t9.*x.*(1.0./2.0),t9.*x.*y.*(3.0./2.0),t4.*t9.*x.*(1.0./2.0)];


function Fxy = Fxy_legendre_3(x,y)
%FXY_LEGENDRE_3
%    FXY = FXY_LEGENDRE_3(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    14-Apr-2015 22:28:40

t2 = y.^2;
t3 = t2.*(1.5e1./2.0);
t4 = t3-3.0./2.0;
t5 = x.^2;
t6 = t5.*(1.5e1./2.0);
t7 = t6-3.0./2.0;
Fxy = [0.0,0.0,0.0,0.0,0.0,1.0,y.*3.0,t4,0.0,x.*3.0,x.*y.*9.0,t4.*x.*3.0,0.0,t7,t7.*y.*3.0,t4.*t7];


function Fxx = Fxx_legendre_3(x,y)
%FXX_LEGENDRE_3
%    FXX = FXX_LEGENDRE_3(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    14-Apr-2015 22:28:40

t2 = y.^2;
t3 = t2.*5.0;
t4 = t3-3.0;
Fxx = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,y.*3.0,t2.*(9.0./2.0)-3.0./2.0,t4.*y.*(3.0./2.0),x.*1.5e1,x.*y.*1.5e1,x.*(t2.*(3.0./2.0)-1.0./2.0).*1.5e1,t4.*x.*y.*(1.5e1./2.0)];


function Fyy = Fyy_legendre_3(x,y)
%FYY_LEGENDRE_3
%    FYY = FYY_LEGENDRE_3(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    14-Apr-2015 22:28:40

t2 = x.^2;
t3 = t2.*5.0;
t4 = t3-3.0;
Fyy = [0.0,0.0,3.0,y.*1.5e1,0.0,0.0,x.*3.0,x.*y.*1.5e1,0.0,0.0,t2.*(9.0./2.0)-3.0./2.0,y.*(t2.*(3.0./2.0)-1.0./2.0).*1.5e1,0.0,0.0,t4.*x.*(3.0./2.0),t4.*x.*y.*(1.5e1./2.0)];
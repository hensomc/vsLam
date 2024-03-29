function [F,Fx,Fy,Fxy,Fxx,Fyy] = FdF_legendre5(nDeg,x,y)

F   = F_legendre_5(x,y);
Fx  = Fx_legendre_5(x,y);
Fy  = Fy_legendre_5(x,y);
Fxy = Fxy_legendre_5(x,y);
Fxx = Fxx_legendre_5(x,y);
Fyy = Fyy_legendre_5(x,y);

function F = F_legendre_5(x,y)
%F_LEGENDRE_5
%    F = F_LEGENDRE_5(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:30:36

t2 = y.^2;
t3 = t2.^2;
t4 = t2.*(3.0./2.0);
t5 = t4-1.0./2.0;
t6 = t2.*5.0;
t7 = t6-3.0;
t8 = t3.*(3.5e1./8.0);
t15 = t2.*(1.5e1./4.0);
t9 = t8-t15+3.0./8.0;
t10 = t3.*6.3e1;
t16 = t2.*7.0e1;
t11 = t10-t16+1.5e1;
t12 = x.^2;
t13 = t12.*(3.0./2.0);
t14 = t13-1.0./2.0;
t17 = t12.*5.0;
t18 = t17-3.0;
t19 = t12.^2;
t20 = t19.*(3.5e1./8.0);
t22 = t12.*(1.5e1./4.0);
t21 = t20-t22+3.0./8.0;
t23 = t19.*6.3e1;
t25 = t12.*7.0e1;
t24 = t23-t25+1.5e1;
F = [1.0,y,t5,t7.*y.*(1.0./2.0),t9,t11.*y.*(1.0./8.0),x,x.*y,t5.*x,t7.*x.*y.*(1.0./2.0),t9.*x,t11.*x.*y.*(1.0./8.0),t14,t14.*y,t5.*t14,t7.*t14.*y.*(1.0./2.0),t9.*t14,t11.*t14.*y.*(1.0./8.0),t18.*x.*(1.0./2.0),t18.*x.*y.*(1.0./2.0),t5.*t18.*x.*(1.0./2.0),t7.*t18.*x.*y.*(1.0./4.0),t9.*t18.*x.*(1.0./2.0),t11.*t18.*x.*y.*(1.0./1.6e1),t21,t21.*y,t5.*t21,t7.*t21.*y.*(1.0./2.0),t9.*t21,t11.*t21.*y.*(1.0./8.0),t24.*x.*(1.0./8.0),t24.*x.*y.*(1.0./8.0),t5.*t24.*x.*(1.0./8.0),t7.*t24.*x.*y.*(1.0./1.6e1),t9.*t24.*x.*(1.0./8.0),t11.*t24.*x.*y.*(1.0./6.4e1)];

function Fx = Fx_legendre_5(x,y)
%FX_LEGENDRE_5
%    FX = FX_LEGENDRE_5(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:30:36

t2 = y.^2;
t3 = t2.^2;
t4 = t2.*(3.0./2.0);
t5 = t4-1.0./2.0;
t6 = t2.*5.0;
t7 = t6-3.0;
t8 = t3.*(3.5e1./8.0);
t15 = t2.*(1.5e1./4.0);
t9 = t8-t15+3.0./8.0;
t10 = t3.*6.3e1;
t16 = t2.*7.0e1;
t11 = t10-t16+1.5e1;
t12 = x.^2;
t13 = t12.*(1.5e1./2.0);
t14 = t13-3.0./2.0;
t17 = t12.*x.*(3.5e1./2.0);
t18 = x.*1.4e2;
t22 = t12.*x.*2.52e2;
t19 = t18-t22;
t20 = t12.^2;
t21 = t20.*(6.3e1./8.0);
t23 = t19.*x.*(1.0./8.0);
t24 = t12.*(3.5e1./4.0);
Fx = [0.0,0.0,0.0,0.0,0.0,0.0,1.0,y,t5,t7.*y.*(1.0./2.0),t9,t11.*y.*(1.0./8.0),x.*3.0,x.*y.*3.0,t5.*x.*3.0,t7.*x.*y.*(3.0./2.0),t9.*x.*3.0,t11.*x.*y.*(3.0./8.0),t14,t14.*y,t5.*t14,t7.*t14.*y.*(1.0./2.0),t9.*t14,t11.*t14.*y.*(1.0./8.0),t17-x.*(1.5e1./2.0),y.*(t17-x.*(1.5e1./2.0)),t5.*(t17-x.*(1.5e1./2.0)),t7.*y.*(t17-x.*(1.5e1./2.0)).*(1.0./2.0),t9.*(t17-x.*(1.5e1./2.0)),t11.*y.*(t17-x.*(1.5e1./2.0)).*(1.0./8.0),t12.*(-3.5e1./4.0)+t21-t19.*x.*(1.0./8.0)+1.5e1./8.0,y.*(t21-t23-t24+1.5e1./8.0),t5.*(t21-t23-t24+1.5e1./8.0),t7.*y.*(t21-t23-t24+1.5e1./8.0).*(1.0./2.0),t9.*(t21-t23-t24+1.5e1./8.0),t11.*y.*(t21-t23-t24+1.5e1./8.0).*(1.0./8.0)];

function Fy = Fy_legendre_5(x,y)
%FY_LEGENDRE_5
%    FY = FY_LEGENDRE_5(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:30:36

t2 = y.^2;
t3 = t2.*(1.5e1./2.0);
t4 = t3-3.0./2.0;
t5 = t2.*y.*(3.5e1./2.0);
t6 = y.*1.4e2;
t13 = t2.*y.*2.52e2;
t7 = t6-t13;
t8 = t2.^2;
t9 = t8.*(6.3e1./8.0);
t10 = x.^2;
t11 = t10.*(3.0./2.0);
t12 = t11-1.0./2.0;
t14 = t7.*y.*(1.0./8.0);
t15 = t2.*(3.5e1./4.0);
t16 = t10.*5.0;
t17 = t16-3.0;
t18 = t10.^2;
t19 = t18.*(3.5e1./8.0);
t21 = t10.*(1.5e1./4.0);
t20 = t19-t21+3.0./8.0;
t22 = t18.*6.3e1;
t24 = t10.*7.0e1;
t23 = t22-t24+1.5e1;
Fy = [0.0,1.0,y.*3.0,t4,t5-y.*(1.5e1./2.0),t2.*(-3.5e1./4.0)+t9-t7.*y.*(1.0./8.0)+1.5e1./8.0,0.0,x,x.*y.*3.0,t4.*x,x.*(t5-y.*(1.5e1./2.0)),x.*(t9-t14-t15+1.5e1./8.0),0.0,t12,t12.*y.*3.0,t4.*t12,t12.*(t5-y.*(1.5e1./2.0)),t12.*(t9-t14-t15+1.5e1./8.0),0.0,t17.*x.*(1.0./2.0),t17.*x.*y.*(3.0./2.0),t4.*t17.*x.*(1.0./2.0),t17.*x.*(t5-y.*(1.5e1./2.0)).*(1.0./2.0),t17.*x.*(t9-t14-t15+1.5e1./8.0).*(1.0./2.0),0.0,t20,t20.*y.*3.0,t4.*t20,t20.*(t5-y.*(1.5e1./2.0)),t20.*(t9-t14-t15+1.5e1./8.0),0.0,t23.*x.*(1.0./8.0),t23.*x.*y.*(3.0./8.0),t4.*t23.*x.*(1.0./8.0),t23.*x.*(t5-y.*(1.5e1./2.0)).*(1.0./8.0),t23.*x.*(t9-t14-t15+1.5e1./8.0).*(1.0./8.0)];


function Fxy = Fxy_legendre_5(x,y)
%FXY_LEGENDRE_5
%    FXY = FXY_LEGENDRE_5(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:30:37

t2 = y.^2;
t3 = t2.*(1.5e1./2.0);
t4 = t3-3.0./2.0;
t5 = t2.*y.*(3.5e1./2.0);
t6 = y.*1.4e2;
t13 = t2.*y.*2.52e2;
t7 = t6-t13;
t8 = t2.^2;
t9 = t8.*(6.3e1./8.0);
t10 = x.^2;
t11 = t10.*(1.5e1./2.0);
t12 = t11-3.0./2.0;
t14 = t7.*y.*(1.0./8.0);
t15 = t2.*(3.5e1./4.0);
t16 = t10.*x.*(3.5e1./2.0);
t17 = x.*1.4e2;
t21 = t10.*x.*2.52e2;
t18 = t17-t21;
t19 = t10.^2;
t20 = t19.*(6.3e1./8.0);
t22 = t18.*x.*(1.0./8.0);
t23 = t10.*(3.5e1./4.0);
Fxy = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,y.*3.0,t4,t5-y.*(1.5e1./2.0),t2.*(-3.5e1./4.0)+t9-t7.*y.*(1.0./8.0)+1.5e1./8.0,0.0,x.*3.0,x.*y.*9.0,t4.*x.*3.0,x.*(t5-y.*(1.5e1./2.0)).*3.0,x.*(t9-t14-t15+1.5e1./8.0).*3.0,0.0,t12,t12.*y.*3.0,t4.*t12,t12.*(t5-y.*(1.5e1./2.0)),t12.*(t9-t14-t15+1.5e1./8.0),0.0,t16-x.*(1.5e1./2.0),y.*(t16-x.*(1.5e1./2.0)).*3.0,t4.*(t16-x.*(1.5e1./2.0)),(t16-x.*(1.5e1./2.0)).*(t5-y.*(1.5e1./2.0)),(t16-x.*(1.5e1./2.0)).*(t9-t14-t15+1.5e1./8.0),0.0,t10.*(-3.5e1./4.0)+t20-t18.*x.*(1.0./8.0)+1.5e1./8.0,y.*(t20-t22-t23+1.5e1./8.0).*3.0,t4.*(t20-t22-t23+1.5e1./8.0),(t5-y.*(1.5e1./2.0)).*(t20-t22-t23+1.5e1./8.0),(t9-t14-t15+1.5e1./8.0).*(t20-t22-t23+1.5e1./8.0)];


function Fxx = Fxx_legendre_5(x,y)
%FXX_LEGENDRE_5
%    FXX = FXX_LEGENDRE_5(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:30:37

t2 = y.^2;
t3 = t2.^2;
t4 = t2.*5.0;
t5 = t4-3.0;
t6 = t3.*6.3e1;
t15 = t2.*7.0e1;
t7 = t6-t15+1.5e1;
t8 = x.^2;
t9 = t8.*(1.05e2./2.0);
t10 = t9-1.5e1./2.0;
t11 = t2.*(3.0./2.0);
t12 = t11-1.0./2.0;
t13 = t3.*(3.5e1./8.0);
t22 = t2.*(1.5e1./4.0);
t14 = t13-t22+3.0./8.0;
t16 = t8.*7.56e2;
t17 = t16-1.4e2;
t18 = t17.*x.*(1.0./8.0);
t19 = t8.*x.*6.3e1;
t21 = x.*3.5e1;
t20 = t18+t19-t21;
Fxx = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,y.*3.0,t2.*(9.0./2.0)-3.0./2.0,t5.*y.*(3.0./2.0),t2.*(-4.5e1./4.0)+t3.*(1.05e2./8.0)+9.0./8.0,t7.*y.*(3.0./8.0),x.*1.5e1,x.*y.*1.5e1,t12.*x.*1.5e1,t5.*x.*y.*(1.5e1./2.0),t14.*x.*1.5e1,t7.*x.*y.*(1.5e1./8.0),t10,t10.*y,t10.*t12,t5.*t10.*y.*(1.0./2.0),t10.*t14,t7.*t10.*y.*(1.0./8.0),t20,t20.*y,t12.*t20,t5.*t20.*y.*(1.0./2.0),t14.*t20,t7.*t20.*y.*(1.0./8.0)];


function Fyy = Fyy_legendre_5(x,y)
%FYY_LEGENDRE_5
%    FYY = FYY_LEGENDRE_5(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:30:38

t2 = y.^2;
t3 = t2.*(1.05e2./2.0);
t4 = t3-1.5e1./2.0;
t5 = t2.*7.56e2;
t6 = t5-1.4e2;
t7 = t6.*y.*(1.0./8.0);
t8 = t2.*y.*6.3e1;
t13 = y.*3.5e1;
t9 = t7+t8-t13;
t10 = x.^2;
t11 = t10.*(3.0./2.0);
t12 = t11-1.0./2.0;
t14 = t10.*5.0;
t15 = t14-3.0;
t16 = t10.^2;
t17 = t16.*(3.5e1./8.0);
t19 = t10.*(1.5e1./4.0);
t18 = t17-t19+3.0./8.0;
t20 = t16.*6.3e1;
t22 = t10.*7.0e1;
t21 = t20-t22+1.5e1;
Fyy = [0.0,0.0,3.0,y.*1.5e1,t4,t9,0.0,0.0,x.*3.0,x.*y.*1.5e1,t4.*x,t9.*x,0.0,0.0,t10.*(9.0./2.0)-3.0./2.0,t12.*y.*1.5e1,t4.*t12,t9.*t12,0.0,0.0,t15.*x.*(3.0./2.0),t15.*x.*y.*(1.5e1./2.0),t4.*t15.*x.*(1.0./2.0),t9.*t15.*x.*(1.0./2.0),0.0,0.0,t10.*(-4.5e1./4.0)+t16.*(1.05e2./8.0)+9.0./8.0,t18.*y.*1.5e1,t4.*t18,t9.*t18,0.0,0.0,t21.*x.*(3.0./8.0),t21.*x.*y.*(1.5e1./8.0),t4.*t21.*x.*(1.0./8.0),t9.*t21.*x.*(1.0./8.0)];

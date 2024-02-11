function [F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_Sym(ipoltype,Nmax)
%% Purpose:  Generate symbolic F and its derivatives for 2 polynomial types
%  ipoltype: (1=power series, 2=legendre)

if(ipoltype == 1)
    [F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_power_sym(Nmax);
end

if(ipoltype == 2)
    [F,dFx,dFy,dFxy,dFxx,dFyy]=FdF_legendre_sym(Nmax);
end

end
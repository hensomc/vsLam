function [ F dFx dFy dFxy dFxx dFyy ] = FdF_legendre_sym( nDeg )

syms x y s;
% [ Px, dPx, dPxx ] = legendre_poly_sym ( nDeg )
% [ Py, dPy, dPyy ] = legendre_poly_sym ( nDeg )
% Py=subs(Py, x, y); dPy=subs(dPy, x, y); dPyy=subs(dPyy, x, y);

[Ps,Px]=Legendre_basis_sym_x(nDeg);

Px=Ps;
Px=subs(Px, s, x);

Py=Ps;
Py=subs(Py, s, y);


dPx=diff(Px,x,1);
dPxx=diff(Px,x,2);

dPy=diff(Py,y,1);
dPyy=diff(Py,y,2);

F    = kron( Px, Py);
dFx  = kron(dPx, Py);
dFy  = kron( Px,dPy);
dFxy = kron(dPx,dPy);
dFxx = kron(dPxx,Py);
dFyy = kron( Px,dPyy);

end


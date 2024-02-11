function FdF_legendre_write_sym(nDeg)
%% creates an m-file function for 2-D legendre polynomial of degree nDeg

[ F Fx Fy Fxy Fxx Fyy ] = FdF_legendre_sym( nDeg );

% FdF.F=F;
% FdF.Fx=Fx;
% FdF.Fy=Fy;
% FdF.Fxy=Fxy;
% FdF.Fxx=Fxx;
% FdF.Fyy=Fyy;

fname = strcat('F_legendre_',int2str(nDeg),'.m');
matlabFunction(F,'file',fname)

fname = strcat('Fx_legendre_',int2str(nDeg),'.m');
matlabFunction(Fx,'file',fname)

fname = strcat('Fy_legendre_',int2str(nDeg),'.m');
matlabFunction(Fy,'file',fname)

fname = strcat('Fxy_legendre_',int2str(nDeg),'.m');
matlabFunction(Fxy,'file',fname)

fname = strcat('Fxx_legendre_',int2str(nDeg),'.m');
matlabFunction(Fxx,'file',fname)

fname = strcat('Fyy_legendre_',int2str(nDeg),'.m');
matlabFunction(Fyy,'file',fname)


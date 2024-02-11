function FdF_genSymLengendre()
% Notes
% ipoltype = 2 = legendre
% p = 8 = polynomial degree
% leave 'file' arg as-is
% change 'F_legendre_x' to 'F_legendre_p
% modify system command to copy single files into one



syms x y

[F,Fx,Fy,Fxy,Fxx,Fyy]=FdF_Sym(2,11)
matlabFunction(F,'file','F_legendre_11','vars',[x,y])
matlabFunction(Fx,'file','Fx_legendre_11','vars',[x,y])
matlabFunction(Fy,'file','Fy_legendre_11','vars',[x,y])
matlabFunction(Fxy,'file','Fxy_legendre_11','vars',[x,y])
matlabFunction(Fxx,'file','Fxx_legendre_11','vars',[x,y])
matlabFunction(Fyy,'file','Fyy_legendre_11','vars',[x,y])

%fails - system('copy *_legendre_8 FdF_legendre7.m');



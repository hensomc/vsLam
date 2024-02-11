function [B] = getRitzBasisVectors(ipoltype,pDeg,psiI,etaI,plot_flag)

if( plot_flag > 0 )
    disp('Calculating basis functions...');
end

% Evaluate polynomial basis functions for integration points psiI and etaI
[F,Fx,Fy,Fxy,Fxx,Fyy]=FdF(ipoltype,pDeg,psiI,etaI);

% Displacement basis functions
B.U=F;
B.Ux=Fx;
B.Uy=Fy;
B.Uxy=Fxy;
B.Uxx=Fxx;
B.Uyy=Fyy;

B.V=F;
B.Vx=Fx;
B.Vy=Fy;
B.Vxy=Fxy;
B.Vxx=Fxx;
B.Vyy=Fyy;

B.W=F;
B.Wx=Fx;
B.Wy=Fy;
B.Wxy=Fxy;
B.Wxx=Fxx;
B.Wyy=Fyy;

B.P=F;
B.Px=Fx;
B.Py=Fy;
B.Pxy=Fxy;
B.Pxx=Fxx;
B.Pyy=Fyy;

B.R=F;
B.Rx=Fx;
B.Ry=Fy;
B.Rxy=Fxy;
B.Rxx=Fxx;
B.Ryy=Fyy;



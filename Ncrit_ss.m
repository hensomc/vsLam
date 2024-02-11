function Ncrit=Ncrit_ss(D,m,n,a,b,Nx,Ny)
%Ncrit_ss - calculates critical load for orthotropic plate loaded in
%biaxial compression
%
% Reference: Whitney pp. 104
%
% Inputs:
% D         = D-matrix plate bending rigidities
% m,n       = # buckling waves in x-dir, y-dir
% a,b       = plate dimensions
% Nx,Ny     = in-plane loads
% 
% Outputs:
% Ncrit     = critical load

% load ratio
k=Ny/Nx;

% plate aspect ratio
R=a/b;

Ncrit = (pi^2)*( D(1,1)*(m^4) + 2.*( D(1,2) + 2.*D(3,3) )*(m^2)*(n^2)*(R^2) + D(2,2)*(n^4)*(R^4) )/ ...
               (a^2*( m^2 + k*(n^2)*(R^2) ));
         
end
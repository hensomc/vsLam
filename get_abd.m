function abd = get_abd(panXY,psi,eta,lamDATA)

aij=[0 0 0; 0 0 0; 0 0 0];
bij=[0 0 0; 0 0 0; 0 0 0];
dij=[0 0 0; 0 0 0; 0 0 0];

% Recover panel geometry
% a=panXY(2,1) - panXY(1,1);
% b=panXY(3,2) - panXY(2,2);
a=2; b=2;  % use isoparametric panel geom to get ply angles is s-t space
    
% Recover laminate data
imat=lamDATA.imat;
theta0=lamDATA.theta0;
theta1=lamDATA.theta1;
thk=lamDATA.thk;
mat_names=lamDATA.mat_names;
phi_rot=lamDATA.phi_rot;

e1   = lamDATA.e1;
e2   = lamDATA.e2;
nu12 = lamDATA.nu12;
g12  = lamDATA.g12;
rho  = lamDATA.rho;

% Layer thickness functions
if(lamDATA.tDeg > 0)
  thk=thk_legendre(lamDATA.tDeg,psi,eta)*lamDATA.thk_coeff';
end
%thk=layer_thk*thk;
%thk=sum(thk_legendre(1,psi,eta));

% Ply offsets
nply = length(theta0);
t=sum(thk);
z(1) = -t/2.;
for(k=1:nply)
  z(k+1) = z(k)+thk(k);
end


% ABD calcs
for(k=1:nply)
  %qij = get_qij(theta0(k),e1,e2,nu12,g12);
  m=imat(k);
  qij = get_qij(get_theta(a,b,theta0(k),theta1(k),phi_rot,psi,eta,k),e1(m),e2(m),nu12(m),g12(m));
  aij = aij + qij*(z(k+1) - z(k));
  bij = bij + qij*(z(k+1)*z(k+1) - z(k)*z(k));
  dij = dij + qij*(z(k+1)*z(k+1)*z(k+1)-z(k)*z(k)*z(k));
  %disp('theta,psi =');get_theta(a,b,theta0(k),theta1(k),phi_rot,psi,eta,k), psi, eta
end


bij=0.5*bij;
dij=dij/3.;
abd = [aij bij;
       bij dij;];
end

function [eqx, eqy, eqxn, eqyn] = get_equiv_lam_props(panXY, lamDATA)

a=2; b=2;  % use isoparametric panel geom to get ply angles is s-t space

% Lamina properties - assumed const over calcs below
e1   = lamDATA.e1;
e2   = lamDATA.e2;
nu12 = lamDATA.nu12;
g12  = lamDATA.g12;
rho  = lamDATA.rho;


% Properties as f(x)
x=linspace(-a/2, a/2);
y=zeros(length(x));
dx = abs(x(2)-x(1));
eqx(1:4)=0;
for(i=1:length(x))
  abd = get_abd(panXY, x(i), y(i), lamDATA);
  eij = lam_engr_constants(abd, lamDATA.thk);
  ex(i) = eij(1);
  ey(i) = eij(2);
  nuxy(i) = eij(3);
  gxy(i) = eij(4);
  eqx(1)=eqx(1)+ex(i)*dx;
  eqx(2)=eqx(2)+ey(i)*dx;
  eqx(3)=eqx(3)+gxy(i)*dx;
  eqx(4)=eqx(4)+nuxy(i)*dx;
end
for i=1:length(eqx)
    eqx(i) = eqx(i)/a;
end
% Normalized x-dir effective stiffnesses w.r.t. lamina props
eqxn(1) = eqx(1)/e1;
eqxn(2) = eqx(2)/e1;
eqxn(3) = eqx(3)/g12;
eqxn(4) = eqx(4)/nu12;

% Properties as f(y)
y=linspace(-b/2, b/2);
x=zeros(length(y));
dy = abs(y(2)-y(1));
eqy(1:4)=0;
for(i=1:length(y))
  abd = get_abd(panXY, x(i), y(i), lamDATA);
  eij = lam_engr_constants(abd, lamDATA.thk);
  ex(i) = eij(1);
  ey(i) = eij(2);
  nuxy(i) = eij(3);
  gxy(i) = eij(4);
  eqy(1)=eqy(1)+ex(i)*dy;
  eqy(2)=eqy(2)+ey(i)*dy;
  eqy(3)=eqy(3)+gxy(i)*dy;
  eqy(4)=eqy(4)+nuxy(i)*dy;
end
for i=1:length(eqy)
    eqy(i) = eqy(i)/b;
end
% Normalized y-dir effective stiffnesses w.r.t. lamina props
eqyn(1) = eqy(1)/e1;
eqyn(2) = eqy(2)/e1;
eqyn(3) = eqy(3)/g12;
eqyn(4) = eqy(4)/nu12;
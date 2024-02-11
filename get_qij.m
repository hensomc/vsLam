function qij = get_qij(theta,e1,e2,nu12,g12)
%get_qij(theta,e1,e2,nu12,g12) returns 3x3 qij 
%    theta = ply angle
%    e1, e2, nu12, g12 = lamina ply properites

qij=zeros(3,3);

theta = theta * (pi/180.);

nu21 = nu12*(e2/e1);

q11 = e1/(1-nu12*nu21);
q12 = (nu21*e1)/(1-nu12*nu21);
q22 = e2/(1-nu12*nu21);
q66 = g12;

c2=cos(2.*theta);
c4=cos(4.*theta);
s2=sin(2.*theta);
s4=sin(4.*theta);

u1 = (3.*q11 + 3.*q22 + 2.*q12 + 4.*q66)/8.;
u2 = (q11 - q22)/2.;
u3 = (q11 + q22 - 2.*q12 - 4.*q66)/8.;
u4 = (q11 + q22 + 6.*q12 - 4.*q66)/8.;
u5 = (q11 + q22 - 2.*q12 + 4.*q66)/8.;

% Consolidate for speed (see in-line matrix below)
% q11b = u1 + u2*c2 + u3*c4;
% q12b = u4 - u3*c4;
% q16b = 0.5*u2*s2 + u3*s4;
% q22b = u1 - u2*c2 + u3*c4;
% q26b = 0.5*u2*s2 - u3*s4;
% q66b = u5 - u3*c4;

% qij = [q11b q12b q16b; 
%        q12b q22b q26b; 
%        q16b q26b q66b];

qij(1,1) = u1 + u2*c2 + u3*c4;
qij(1,2) = u4 - u3*c4;
qij(1,3) = 0.5*u2*s2 + u3*s4;
qij(2,2) = u1 - u2*c2 + u3*c4;
qij(2,3) = 0.5*u2*s2 - u3*s4;
qij(3,3) = u5 - u3*c4;
qij(2,1) = qij(1,2);
qij(3,1) = qij(1,3);
qij(3,2) = qij(2,3);

end

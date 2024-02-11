function lam_eij = lam_engr_constants(abd, thk);
%lam_engr_constants(abd);
%    abd = laminate abd stiffness matrix
%    thk = ply thickness array

a = abd(1:3,1:3);
b = abd(1:3,4:6);
d = abd(4:6,4:6);

t = sum(thk);
ap = inv(a - b*inv(d)*b);

ex   = 1/(t*ap(1,1));       %Ex
ey   = 1/(t*ap(2,2));       %Ey
nuxy = -ap(1,2)/ap(1,1);    %nuxy
gxy  = 1/(t*ap(3,3));       %gxy

lam_eij = [ex ey nuxy gxy];

end

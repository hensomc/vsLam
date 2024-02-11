function fibK = getFibCurvatureSkew(b,T0,T1,x,phi)
% Returns instantaneous fiber path curvature

T0i= T0;
T0 = T0 * (pi/180.);
T1 = T1 * (pi/180.);

fibK = (2*(T1-T0)/(b*cosd(phi)))*cos( T0 + 2*(T1-T0)*x/(b*cosd(phi)));

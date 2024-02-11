function fibK = getFibCurvature(a,T0,T1,x)
% Returns instantaneous fiber path curvature

T0i= T0;
T0 = T0 * (pi/180.);
T1 = T1 * (pi/180.);

fibK = (2/a)*(T1-T0)*cosd( T0 + (2/a)*(T1-T0)*x );

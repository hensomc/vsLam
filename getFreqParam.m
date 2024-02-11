function [freqParam]=getFreqParam(EgN,panXY,lamDATA,method)
%% Frequency parameter of 1st mode
a=panXY(2,1) - panXY(1,1);

% Laminate thickness
h=sum(lamDATA.thk);

% Laminate stiffnesses
abd=get_abd(panXY, 0, 0, lamDATA);  % need to transform to psi,eta
d=abd(4:6, 4:6);

% Recover lamina properties - use first ply as representative of laminate
e1   = lamDATA.e1(1);
e2   = lamDATA.e2(1);
nu12 = lamDATA.nu12(1);
g12  = lamDATA.g12(1);
rho  = lamDATA.rho(1);

% switch on method calculation
switch method
    
    case 1  % Houmat Method
        freqParam=(1/(2*pi))*(EgN(1)*(a))*sqrt( (mean(lamDATA.rho_ply))/((12.*32.2)*e2) );
        
    case 2  % Timarci Definition
        %freqParam=(EgN(1)*(a^2))*sqrt( (lamDATA.rho_ply*lamDATA.thk')/((12.*32.2)/d(2,2)) );
        Do=e2*h^3;
        freqParam=(1/(2*pi))*( EgN(1)*a^2)*sqrt( (lamDATA.rho_ply*lamDATA.thk')/((12.*32.2)*Do) );
        
end
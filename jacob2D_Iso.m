function J = jacob2D_Iso( psi, eta, panXY )
%JACOB2D_ISO Compute jacobian in 2-D in psi-eta system

%   Detailed explanation goes here
    dx_dpsi = (1/4)*((1+eta)*(panXY(3,1)-panXY(4,1)) + (1-eta)*(panXY(2,1)-panXY(1,1)));
    dx_deta = (1/4)*((1+psi)*(panXY(3,1)-panXY(2,1)) + (1-psi)*(panXY(4,1)-panXY(1,1)));
    dy_dpsi = (1/4)*((1+eta)*(panXY(3,2)-panXY(4,2)) + (1-eta)*(panXY(2,2)-panXY(1,2)));
    dy_deta = (1/4)*((1+psi)*(panXY(3,2)-panXY(2,2)) + (1-psi)*(panXY(4,2)-panXY(1,2)));
    
%     J = [dx_dpsi dx_deta;
%          dy_dpsi dy_deta];
    J = [dx_dpsi dy_dpsi;
         dx_deta dy_deta];   
end


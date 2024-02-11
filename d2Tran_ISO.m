function T = d2Tran_ISO( psi, eta, panXY )
% d2Tran computes a transformation matrix for 2nd derivatives of a 
% function F(x,t) to F(s,t)

% The matrix is ordered as follows
% [F''(s,t)] = [T]*[F''(x,y)]
%
% Apply as:  [F''(x,y)] = inv[T]*[F''(s,t)]

% [F''(x,y)] = [Fxx; Fxy; Fyy]
% [F''(s,t)] = [Fss  Fst  Ftt]

    dx_dpsi = (1/4)*((1+eta)*(panXY(3,1)-panXY(4,1)) + (1-eta)*(panXY(2,1)-panXY(1,1)));
    dx_deta = (1/4)*((1+psi)*(panXY(3,1)-panXY(2,1)) + (1-psi)*(panXY(4,1)-panXY(1,1)));
    dy_dpsi = (1/4)*((1+eta)*(panXY(3,2)-panXY(4,2)) + (1-eta)*(panXY(2,2)-panXY(1,2)));
    dy_deta = (1/4)*((1+psi)*(panXY(3,2)-panXY(2,2)) + (1-psi)*(panXY(4,2)-panXY(1,2)));
    
    T = [dx_dpsi^2        2*dx_dpsi*dy_dpsi                dy_dpsi^2
         dx_dpsi*dx_deta  dx_deta*dy_dpsi+dx_dpsi*dy_deta  dy_dpsi*dy_deta
         dx_deta^2        2*dx_deta*dy_deta                dy_deta^2];
     
 
%      J = [dx_dpsi dy_dpsi;
%          dx_deta dy_deta];
     J=jacob2D_Iso( psi, eta, panXY );
     
     T=[J(1,1)^2       2*J(1,1)*J(1,2)              J(1,2)^2
        J(2,1)^2       2*J(2,1)*J(2,2)              J(2,2)^2
        J(1,1)*J(2,1)  J(2,1)*J(1,2)+J(1,1)*J(2,2)  J(1,2)*J(2,2)];

end


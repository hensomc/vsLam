function T = d1Tran_ISO( psi, eta, panXY )
% d1Tran computes a transformation matrix for 1st derivatives of a 
% function F(x,t) to F(s,t)

% The matrix is ordered as follows
% [F'(s,t)] = [T]*[F'(x,y)]
%
% Apply as:  [F''(x,y)] = inv[T]*[F''(s,t)]

% [F'(x,y)] = [Fx; Fy]
% [F'(s,t)] = [Fs  Ft]

    dx_dpsi = (1/4)*((1+eta)*(panXY(3,1)-panXY(4,1)) + (1-eta)*(panXY(2,1)-panXY(1,1)));
    dx_deta = (1/4)*((1+psi)*(panXY(3,1)-panXY(2,1)) + (1-psi)*(panXY(4,1)-panXY(1,1)));
    dy_dpsi = (1/4)*((1+eta)*(panXY(3,2)-panXY(4,2)) + (1-eta)*(panXY(2,2)-panXY(1,2)));
    dy_deta = (1/4)*((1+psi)*(panXY(3,2)-panXY(2,2)) + (1-psi)*(panXY(4,2)-panXY(1,2)));
    
    T = [dx_dpsi   dy_dpsi
         dx_deta   dy_deta];

end


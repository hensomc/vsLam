function panWt = panelWeightVthk2( lamDATA,panXY,psiI,etaI,gl_wts,intMethod, plot_flag )

%% calculates panel Weight for variable thickness
%
%  Inputs:
%   panXY - array of panel corner points

layer_thk = getLayerThk(lamDATA.tDeg,psiI,etaI,lamDATA,plot_flag);

lamDATA.rho_ply; % Will need to compute rho_avg(psi,eta) for mult layers

% Integrate thickness over plate surface area to obtain weight
panWt=0;
for i=1:length(psiI);
    gl_wt=1.0;
    if(intMethod == 2)  %gaussian quadrature weights
        gl_wt=gl_wts(i);
    end
    
    % Jacobian
    J=jacob2D_Iso( psiI(i), etaI(i), panXY );
    
    rho=(lamDATA.rho_ply); % need a compute a variable rho value (same as layer_thk - maybe do both at once vRhoThk

    panWt = panWt + layer_thk(i)*det(J)*rho*gl_wt;

end

if( plot_flag > 0 )
    disp('Panel Weight = '); panWt
end
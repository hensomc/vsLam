function panWtVthk = panelWeightVthk(solDATA,lamDATA,panXY )
%% calculates panel Weight for variable thickness
%
%  Inputs:
%   panXY - array of panel corner points

% Simplified weight calc for constant thickness layers
if(lamDATA.tDeg == 0)
    %panWtVthk = panelWeight(lamDATA.thk,lamDATA.rho_ply,panXY);
       
    % Convert to 3D by add 0 z-coordinate value
    panXY3D=[panXY(1,:) 0
        panXY(2,:) 0
        panXY(3,:) 0
        panXY(4,:) 0];
    
    area = 0.5*cross(panXY3D(3,:)-panXY3D(1,:),panXY3D(4,:)-panXY3D(2,:));
    area=area(3);
    
    panWtVthk=(lamDATA.thk*lamDATA.rho_ply')*area;
    
    return
else

    % Recover solution parameters
    [ipoltype,pDeg,Mg,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);
    
    % Define Integration points
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(pDeg,intMethod,panXY,plot_flag);
    
    % Polynomial thickness distribution
    vThk = thk_legendre(lamDATA.tDeg,psiI,etaI)*lamDATA.thk_coeff';
    
    lamDATA.rho_ply;
    
    % Integrate thickness over plate surface area to obtain weight
    panWtVthk=0;
    for i=1:length(psiI);
        gl_wt=1.0;
        if(intMethod == 2)  %gaussian quadrature weights
            gl_wt=gl_wts(i);
        end
        
        % Jacobian
        J=jacob2D_Iso( psiI(i), etaI(i), panXY );
        
        rho=(lamDATA.rho_ply); % need to compute a variable rho value (same as vThk - maybe do both at once vRhoThk
        
        panWtVthk = panWtVthk + vThk(i)*det(J)*rho*gl_wt;
        
    end   
end

if(plot_flag > 0)
    disp('Panel Weight = ');panWtVthk
end
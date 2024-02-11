function [ meshDATA, intDATA ] = meshGen( panXY, solDATA )
%meshGen Generates mesh of points for Ritz plate analyses
%   Detailed explanation goes here

% Recover solution parameters
[ipoltype,pDeg,Mg,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);


[psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(pDeg);
[xBCMesh,yBCMesh] = meshTransform(panXY, psiBCMesh, etaBCMesh);

%% Ritz super-element nodal-mesh
% if(RitzSE==1)
%     [psiDomainSE,etaDomainSE]=bcMeshIsoSE(pDeg);
%     xMeshSE=psiDomainSE;    yMeshSE=etaDomainSE;
% end

%% Define Integration points
[psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(pDeg,intMethod,panXY,plot_flag);


%% Return meshDATA
meshDATA.psiBCMesh = psiBCMesh;
meshDATA.etaBCMesh = etaBCMesh;
meshDATA.isoGridID = isoGridID;

meshDATA.xBCMesh = xBCMesh;
meshDATA.yBCMesh = yBCMesh;

intDATA.psiIMesh = psiIMesh;
intDATA.etaIMesh = etaIMesh;
intDATA.psiI = psiI;
intDATA.etaI = etaI;
% intDATA.xIMesh = xIMesh;
% intDATA.yIMesh = yIMesh;
intDATA.gl_wts = gl_wts;

end
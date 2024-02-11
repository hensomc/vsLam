function[psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(pDeg)
%% defines a mesh of Grids in isoparametric coords {psi}&{eta}
nGridPsi = pDeg+1;
nGridEta = pDeg+1;
psiDomain=linspace(-1,1,nGridPsi);
etaDomain=linspace(-1,1,nGridEta);
[psiBCMesh,etaBCMesh]=meshgrid(psiDomain,etaDomain);
isoGridID=[ (1:nGridPsi*nGridEta)' psiBCMesh(:) etaBCMesh(:) ];

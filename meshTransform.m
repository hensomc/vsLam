function [xMesh,yMesh] = meshTransform(panXY, psiMesh, etaMesh)

% Creates an x-y mesh from psi-eta mesh using bi-linear interpolation
% Reddy FEM text pg 252
% panXY = panel corner pt coords in x-y space
% psiMesh, etaMesh = nxn matrices of psi-eta mesh points, nondimensional
%                    coordinates

dim_psi=size(psiMesh);
dim_eta=size(etaMesh);

psi1=reshape(psiMesh,[1,dim_psi(1)*dim_psi(2)]);
eta1=reshape(etaMesh,[1,dim_eta(1)*dim_eta(2)]);

xmesh1=(1/4)*((1-psi1).*(1-eta1)*panXY(1,1) + (1+psi1).*(1-eta1)*panXY(2,1) + (1+psi1).*(1+eta1)*panXY(3,1) + (1-psi1).*(1+eta1)*panXY(4,1));
ymesh1=(1/4)*((1-psi1).*(1-eta1)*panXY(1,2) + (1+psi1).*(1-eta1)*panXY(2,2) + (1+psi1).*(1+eta1)*panXY(3,2) + (1-psi1).*(1+eta1)*panXY(4,2));

xMesh=reshape(xmesh1,[dim_psi(1),dim_psi(2)]);
yMesh=reshape(ymesh1,[dim_eta(1),dim_eta(2)]);
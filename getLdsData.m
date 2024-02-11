function [ Nx, Ny, Nxy, ptLds, ptLoc, p0 ] = getLdsData( ldsDATA )
%setLdsData Sets the external panel loads data

%  ldsDATA = returned data structure

Nx=ldsDATA.Nx;
Ny=ldsDATA.Ny;
Nxy=ldsDATA.Nxy;

ptLds=ldsDATA.ptLds;
ptLoc=ldsDATA.ptLoc;

p0=ldsDATA.p0;
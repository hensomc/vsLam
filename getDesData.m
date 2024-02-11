function [ XL, XU, X0, Xstart, Rallow, drespType ] = getDesData( desDATA )
%setDesData Sets design obtmization data structure

%  desDATA = returned data structure

XL = desDATA.XL;
XU = desDATA.XU;
X0 = desDATA.X0;
Xstart = desDATA.Xstart;

Rallow = desDATA.Rallow;

drespType = desDATA.drespType;
function [ xc, yc ] = ISO_st_to_xy2(psiI, etaI, panXY )
%ISO_st_to_xy2 Transforms [psiI,etaI] to [xc,yc] using 4 point
%interpolation
%   Detailed explanation goes here

    N(:,1)=(1/4)*(1-psiI).*(1-etaI);
    N(:,2)=(1/4)*(1+psiI).*(1-etaI);
    N(:,3)=(1/4)*(1+psiI).*(1+etaI);
    N(:,4)=(1/4)*(1-psiI).*(1+etaI);
    
    xc=0;yc=0;
    for j=1:4
        xc = xc + N(:,j)*panXY(j,1);
        yc = yc + N(:,j)*panXY(j,2);
    end
    xc = xc';
    yc = yc';
end
function [Kq,Mq,Fq,T] = applyBC( panXY,solDATA,edgeBC,edgeGID,isoGridID,bcType,ldsDATA,Kcww,Mc,Fc,plot_flag )

pDeg=solDATA.pDeg;
ipoltype=solDATA.ipoltype;
isoltype=solDATA.isoltype;
bcMethod=solDATA.bcMethod;

if( bcMethod == 2)
    [T]=nullBC_wpr(ipoltype,isoltype,pDeg,panXY,edgeBC,edgeGID,isoGridID,bcType,ldsDATA,plot_flag);
    Kq=T'*Kcww*T;
    Mq=T'*Mc*T;
    Fq=T'*Fc;   % Point loads
elseif( bcMethod == 1)
    [Kspr]=getKspr2(ipoltype,pDeg,bcType,bcSpring,edgeBC,edgeGID,isoGridID);
    Kq=Kcww + Kspr;
    Mq=Mc;
    Fq=Fc;
    T=0;
end

if( plot_flag>0)
    disp({'size of force vector Fc = ',int2str(size(Fc))});
    disp({'size of generalized force vector Fq = ',int2str(size(Fq))});
end

KM_matrix_qualities(Kq,Mq,plot_flag);

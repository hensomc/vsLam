function [T]=nullBC_edge(ipoltype,pDeg,panXY,edgeBC,edgeGID,isoGridID,bcType)
%% Apply null space to edge boundary conditions
%symbolic evaluation
if(ipoltype == 1)
    syms x y;
    [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF_power_sym(pDeg);
end


%Define BC transformation matrix
IBC=0;
%for ibc=1:length(BC)
for i=1:length(edgeBC)
    
    BC=edgeGID(i,:);
    BCpsieta=[isoGridID(BC,2) isoGridID(BC,3)];
    
    if( bcType(i) ~= 'f' )  % skip if edge is free
        K1=0; K2=0; K3=0;
        for ibc=1:length(BC)
            
            IBC=IBC+1; %  Note: IBC and ibc are 2 variables
            
            psi = BCpsieta(ibc,1);
            eta = BCpsieta(ibc,2);
            
            [xc,yc]=ISO_st_to_xy(psi,eta,panXY);
            if(ipoltype == 2 || ipoltype == 3)
                [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF(ipoltype,pDeg,psi,eta);
                %[W,dWx,dWy,dWxy,dWxx,dWyy]=FdF(ipoltype,pDeg,xc,yc);
            end
            
            % constrain w=0
            if(ipoltype == 1)
                F1(IBC,:)=double(subs(W,[x y],[psi eta]));
                %F1(IBC,:)=double(subs(W,[x y], [xc yc]));
            elseif(ipoltype == 2 || ipoltype == 3)
                F1(IBC,:)=W;
            end
            
            % constrain dw/dx=0
            if( bcType(i) == 'c' )
                IBC=IBC+1;
                if(ipoltype == 1)
                    %F1(IBC,:)=double(subs(dWx,[x y],[xc yc]));
                    F1(IBC,:)=double(subs(dWx,[x y],[psi eta]));
                elseif(ipoltype == 2)
                    F1(IBC,:)=dWx;
                end
            end
            
            % constrain dw/dy=0
            if( bcType(i) == 'c' )
                IBC=IBC+1;
                if(ipoltype == 1)
                    %F1(IBC,:)=double(subs(dWy,[x y],[xc yc]));
                    F1(IBC,:)=double(subs(dWy,[x y],[psi eta]));
                elseif(ipoltype == 2)
                    F1(IBC,:)=dWy;
                end
            end
        end
    end
end
FBC=F1;
%T=null(FBC);
T=FBC;
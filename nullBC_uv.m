function [T]=nullBC_uv(ipoltype,pDeg,panXY,edgeBC,edgeGID,isoGridID,bcType)
%% Apply null space to u-v coordinates for edge boundary conditions

% case: symbolic evaluation
if(ipoltype == 1)
    syms x y;
    [U,dUx,dUy,dUxy,dUxx,dUyy]=FdF_power_sym(pDeg);
    [V,dVx,dVy,dVxy,dVxx,dVyy]=FdF_power_sym(pDeg);
end

zvec=zeros(1,(pDeg+1)^2);
%Define BC transformation matrix
IBC=0;

%for ibc=1:length(BC)
for i=1:length(edgeBC)
    
    BC=edgeGID(i,:);
    BCpsieta=[isoGridID(BC,2) isoGridID(BC,3)];
    
    if( bcType(i) ~= 'f' )  % skip if edge is free
        K1=0; K2=0; K3=0;
        lenBC=length(BC);
        
        % Unit vector for edge
        edgeVec = [BCpsieta(lenBC,1)-BCpsieta(1,1) ...
                   BCpsieta(lenBC,2)-BCpsieta(1,2) ...
                   0];
        uvec = createUnitVector(edgeVec);

        for ibc=1:lenBC;
            
            IBC=IBC+1; %  Note: IBC and ibc are 2 variables
            
            psi = BCpsieta(ibc,1);
            eta = BCpsieta(ibc,2);
            
            [xc,yc]=ISO_st_to_xy(psi,eta,panXY);
            if(ipoltype == 2 || ipoltype == 3)
                [U,dUx,dUy,dUxy,dUxx,dUyy]=FdF(ipoltype,pDeg,psi,eta);
                [V,dVx,dVy,dVxy,dVxx,dVyy]=FdF(ipoltype,pDeg,psi,eta);
            end
            
            % constrain displacement
            if( bcType(i) == 'c' || bcType(i) == 's')
                %disp('RitzBC = ');BC
                IBC=IBC+1;
                % u = 0 for edges parallel to psi (y-dir)
                if( uvec(2) == 1.0 )
                    %disp('u = 0');
                    if(ipoltype == 1)
                        F1(IBC,:)=double(subs(U,[x y],[psi eta]));
                    elseif(ipoltype == 2)
                        F1(IBC,:)=[U zvec];
                        IBC=IBC+1;
                        F1(IBC,:)=[zvec V];
                    end

                % v = 0 for edges parallel to psi (x-dir)
                elseif( uvec(1) == 1.0 )
                    %disp('v = 0');
                    if(ipoltype == 1)
                        F1(IBC,:)=double(subs(V,[x y],[psi eta]));
                    elseif(ipoltype == 2)
                        F1(IBC,:)=[V zvec];
                    end
                end
            end
        end
    end
end
FBC=F1;
T=null(FBC);
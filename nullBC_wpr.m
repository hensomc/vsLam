function [T]=nullBC_wpr(ipoltype,isoltype,pDeg,panXY,edgeBC,edgeGID,isoGridID,bcType,ldsDATA,plot_flag)
%% Apply null space to edge boundary conditions

if(plot_flag > 0)
    fprintf('Generating BC matrix...\n');
end

%symbolic evaluation
if(ipoltype == 1)
    syms x y;
    [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF_power_sym(pDeg);
end

zvec=zeros(1,(pDeg+1)^2);

% recover loads data
[ Nx, Ny, Nxy, ptLds, ptLoc, p0 ] = getLdsData( ldsDATA );


%Define BC transformation matrix
IBC=0;

% Loop over panel edges with defined BCs = [4]
for i=1:length(edgeBC)
    length(edgeBC);
    BC=edgeGID(i,:);
    BCpsieta=[isoGridID(BC,2) isoGridID(BC,3)];
    
    % f = FREE, Skip block if free
    if( bcType(i) ~= 'f' )
        K1=0; K2=0; K3=0;
        lenBC=length(BC);  % lenBC = [nBCpts^2] = (Mg+1)^2
        
        % Unit vector for edge
        edgeVec = [BCpsieta(lenBC,1)-BCpsieta(1,1) ...
                   BCpsieta(lenBC,2)-BCpsieta(1,2) ...
                   0];
        uvec = createUnitVector(edgeVec);

        for ibc=1:lenBC;  % (Mg+1)^2
            
            IBC=IBC+1; %  Note: IBC and ibc are 2 variables
            
            psi = BCpsieta(ibc,1);
            eta = BCpsieta(ibc,2);
            
            [xc,yc]=ISO_st_to_xy(psi,eta,panXY);
            if(ipoltype == 2 || ipoltype == 3)
               [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF(ipoltype,pDeg,psi,eta); %W displacement
                U=W;dUx=dWx;dUy=dWy;dUxy=dWxy;dUxx=dWxx;dUyy=dWyy;  % U displacement
                V=W;dVx=dWx;dVy=dWy;dVxy=dWxy;dVxx=dWxx;dVyy=dWyy;  % V displacement
                P=W;dPx=dWx;dPy=dWy;dPxy=dWxy;dPxx=dWxx;dPyy=dWyy;  % Phi_x rotation
                R=W;dRx=dWx;dRy=dWy;dRxy=dWxy;dRxx=dWxx;dRyy=dWyy;  % Phi_y rotation
                %[W,dWx,dWy,dWxy,dWxx,dWyy]=FdF(ipoltype,pDeg,xc,yc);
            end
            
            % constrain u,v,w=0
            if(ipoltype == 1)
                F1(IBC,:)=double(subs(W,[x y],[psi eta]));
                %F1(IBC,:)=double(subs(W,[x y], [xc yc]));
                
            elseif(ipoltype == 2 || ipoltype == 3)
                
                if(bcType(i) == 'c')
                    
                    if(isoltype==2)  
                        %if(i==1) % hardwired for Nx-only buckling soln
                        if(i==2) % hardwired for NY-only buckling soln
                            F1(IBC,:)= [U zvec zvec zvec zvec];
                            IBC=IBC+1;
                            F1(IBC,:)= [zvec V zvec zvec zvec];
                            IBC=IBC+1;
                            F1(IBC,:)= [zvec zvec W zvec zvec];
                        else
                            IBC=IBC+1;
                            F1(IBC,:)= [zvec zvec W zvec zvec];
                        end
                            
                    %if(isoltype==2 && i==3) % hardwired for Nx-only buckling soln - allows edge to move in load direction
                    %    F1(IBC,:)= [zvec V zvec zvec zvec];
                    %if(isoltype==2 && i==4) % hardwired for Ny-only buckling soln - allows edge to move in load direction
                    %    F1(IBC,:)= [U zvec zvec zvec zvec];
                    %    IBC=IBC+1;
                    %    F1(IBC,:)= [zvec zvec W zvec zvec];
                    else                                        
                        F1(IBC,:)= [U zvec zvec zvec zvec];
                        IBC=IBC+1;
                        F1(IBC,:)= [zvec V zvec zvec zvec];
                        IBC=IBC+1;
                        F1(IBC,:)= [zvec zvec W zvec zvec];
                    end
                    
                elseif(bcType(i) == 's')
                    
                    %if( abs(Nx) > 0)
                        
                    %elseif (abs(Ny) > 0)
                        
                    %if(i==1)  % constrain only 1st edge in-plane: NX Loads
                    if(i==2)  % constrain 2nd edge: NY Loads (TEMP)
                        F1(IBC,:)= [U zvec zvec zvec zvec];
                        IBC=IBC+1;
                        F1(IBC,:)= [zvec V zvec zvec zvec];  
                        IBC=IBC+1;                        
                        F1(IBC,:)= [zvec zvec W zvec zvec];                        
                    else
                        F1(IBC,:)= [zvec zvec W zvec zvec];
                    end
                    
                end
            end
            
            % constrain rotation
            if( bcType(i) == 'c' )
                %disp('RitzBC = ');BC
                IBC=IBC+1;
                
                if(ipoltype == 1)
                    F1(IBC,:)=double(subs(dWx,[x y],[xc yc]));
                    F1(IBC,:)=double(subs(dWx,[x y],[psi eta]));
                elseif(ipoltype == 2)
                    F1(IBC,:)=[zvec zvec zvec P zvec];
                    IBC=IBC+1;
                    F1(IBC,:)=[zvec zvec zvec zvec R];
                end                           
            end
        end
    end
end
FBC=F1;
T=null(FBC);
if(plot_flag > 0)
    disp({'size of FBC=',int2str(size(FBC))})
    disp({'size of T=',int2str(size(T))})
end

function [G]=K_ExtLds(panXY,isoltype,ipoltype,pDeg, intMethod,psiI,etaI,gl_wts,W,W_ptLds,dWx,dWy,isoGridID,dUx,dUy,dVx,dVy,ldsDATA,plot_flag)
%% Geometric stiffness matrix - Potential energy due to external loads (

zvec = zeros(1,(pDeg+1)^2);

% Recover external panel loads
[Nx,Ny,Nxy,ptLds,ptLoc,p0]=getLdsData(ldsDATA);

G=0;

% Panel geometry variables
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);

% In-plane loads
if( abs(Nx)+abs(Ny)+abs(Nxy) > 0 )
    if(isoltype == 2)  %buckling
        if(plot_flag > 0)
            disp('Calculating geometric stiffness matrix...');
        end
        for i=1:length(psiI);
            
            gl_wt=1.0;
            if(intMethod==2)  % gaussian quadrature
                gl_wt=gl_wts(i);
            end
            
            % Transformation for 1st derivatives
            T  = d1Tran_ISO( psiI(i), etaI(i), panXY );
            TI = inv(T);
            
            % Strains due to a bending deflection
            Wx=dWx(i,:);Wy=dWy(i,:);
            %WxWy_sym = Wx'*Wy + Wy'*Wx;
            
            % Transform to psi-eta
            [Wp]=TI*[Wx; Wy];
            Wx=Wp(1,:);
            Wy=Wp(2,:);

            zz=[zvec zvec   zvec   zvec zvec];
            Wx=[zvec zvec dWx(i,:) zvec zvec];
            Wy=[zvec zvec dWy(i,:) zvec zvec];

            
            eps = [zz; zz; zz; zz; Wx; Wy; zz; zz; zz; zz; zz; zz];  % [12 x 5*(pDeg+1)^2]
                                   
            % Evaluate jacobian
            psi=psiI(i); eta=etaI(i);
            J=jacob2D_Iso( psi, eta, panXY );
            JI=inv(J);
            j11=JI(1,1); j12=JI(1,2); j21=JI(2,1); j22=JI(2,2);
            
            T = [0 0 0 0   j11     j12   0 0 0 0 0 0;
                 0 0 0 0   j21     j22   0 0 0 0 0 0;
                 0 0 0 0 j11+j21 j12+j22 0 0 0 0 0 0;
                 0 0 0 0    0       0    0 0 0 0 0 0;
                 0 0 0 0    0       0    0 0 0 0 0 0];  % [5 x 12]
             
            N=[Nx  0  0  0  0;
               0  Ny  0  0  0;
               0   0 Nxy 0  0;
               0   0  0  0  0;
               0   0  0  0  0];  % [5x5]
                
            % PE of inplane loads due to bending deflection
            %G = G + ( Nx*(Wx'*Wx) + Ny*(Wy'*Wy) + 2*Nxy*(Wx'*Wy) )*det(J)*gl_wt;

            G = G + eps'*T'*N'*T*eps*det(J)*gl_wt;

        end
        
        if(isoltype == 2)
            disp({'G_err = ',num2str(max(max(abs(G-G'))))});
        end

    % Integral of panel edge loads {N}*{q} around plate boundaries
    elseif(isoltype == 3)
        if(plot_flag>0)
            disp('Calculating work by external panel loads Nx, Ny, Nxy...');
        end
        ab=[-1  1; -1  1; 1  -1; 1  -1];  % ab integral bounds for edges
        
        % Loop over edges, integrate {N}{d} along each boundary
        for i=1:length(ab)
            [xpt,wt]=lgwt(pDeg+1,ab(i,1),ab(i,2));
            if(i==1)
                ypt(1:pDeg+1)=-1;
                edgeVec=[ab(i,2)-ab(i,1)/2 0 0];
            elseif(i==2)
                ypt=xpt;
                xpt(1:pDeg+1)=1;
                edgeVec=[0 ab(i,2)-ab(i,1)/2 0];
            elseif(i==3)
                ypt(1:pDeg+1)=1;
                edgeVec=[ab(i,2)-ab(i,1)/2 0 0];               
            elseif(i==4)
                ypt=xpt;
                xpt(1:pDeg+1)=-1;
                edgeVec=[0 ab(i,2)-ab(i,1)/2 0];
            end
            edgeVec = createUnitVector(edgeVec);
            
            % Evaluate displacement functions at integ. pts on edge
            [u]=FdF(ipoltype,pDeg,xpt,ypt);
            [v]=[u];

            % Perform integration
            for j=1:length(wt)
                
                 urow = [u(j,:)  zvec];
                 vrow = [zvec  v(j,:)];
                
                % only edge at psi = 1 considered for now
                if(i==2)
                    J=jacob2D_Iso( xpt(j), ypt(j), panXY );
                    G = G + ( Nx*urow+Nxy*vrow )*J(2,2)*wt(j);
                end
            end
        end
        % Needs to be generalized - must be adjusted for Nx or Nxy
        %G = [G zvec];
        %G = [zvec G];
        
        G=G';
    end
end

% Distributed transverse loads
if( abs(p0) > 0.0 )
    if(plot_flag > 0)
        disp('Calculating work by external pressure loads...');
    end
    for i=1:length(psiI);
        gl_wt=1.0;
        if(intMethod==2)  % gaussian quadrature
            gl_wt=gl_wts(i);
        end
        
        psi=psiI(i); eta=etaI(i);
        J=jacob2D_Iso( psi, eta, panXY );
        
        if(isoltype == 3)
            %G = G + (p0 * W)*det(J)*gl_wt;  % Eq(33)
            G= G + [zvec zvec (p0*W(i,:))*det(J)*gl_wt zvec zvec];
        end
    end
    G=G';
%     G_orig=G;
%     S=sum(G);
%     G=S';
%     G=sum(G,1);
%     G=G'
end

% Bending due to transverse concentrated loads
%disp('ptLds = '); ptLds
if( max(any(abs(ptLds))) > 0 )
    if(plot_flag > 0)
        disp('Calculating work by external point loads...');
    end
    num_ptLds = size(ptLds,1);
    for i=1:num_ptLds
        
        [psi,eta]=ISO_XY_to_st_Num(ptLoc(i,1), ptLoc(i,2), panXY);
               
        %G = G + W_ptLds(i,:)'*ptLds(i,3);  % Eq(33)
        G= G + [zvec zvec W_ptLds(i,:)*ptLds(i,3) zvec zvec]';
    end
end

if( plot_flag > 0 )
    disp({'Size of G = ',int2str(size(G))});
end

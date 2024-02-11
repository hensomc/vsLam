function eigB = getEigOrthoRect(bcType,ilam,a,b,abd)
%geEigOrtho()
%    calculates eigenvalaues for orthotropic rectangular plates
%    eigB  = returned Baseline eigenvalues

%    See Whitney, pp 122-130, also Blevins pp 262, 267-268


% plate aspect ratio
R=a/b;
[imat, theta0, theta1, thk, mat_names, phi_rot] = readLam(ilam);
tlam=sum(thk);

% laminate properties

D = abd(4:6,4:6);


% Material properties
[e1,e2,nu12,g12,rho,mtype] = get_mat_props(mat_names);
trho=tlam*rho/(32.2*12.0);

switch bcType
    case {'scsf'}
        %[Ex Ey nuxy Gxy] = lam_engr_constants(abd, thk);
        lam_eij = lam_engr_constants(abd, thk);
        Ex=lam_eij(1); Ey=lam_eij(2); nuxy=lam_eij(3); Gxy = lam_eij(4);
        nu_x = nuxy
        nu_y = nu_x*Ey/Ex;
        Dx=Ex*tlam^3/(12*(1-nu_x*nu_y))
        Dy=Ey*tlam^3/(12*(1-nu_x*nu_y))
        Dk=Gxy*tlam^3/12
        Dxy=Dx*nu_y+2*Dk
        
%         Dx = D(1,1);
%         Dy = D(2,2);
%         Dxy = D(3,3);
%         Dk = D(1,2);
        
        G = [0.597 1.494 2.500];
        H = [-0.0870 1.347 4.658];
        J = [0.471 3.284 7.842];
        for(i=1:3)
            m=i;
            G1 = G(i);
            H1 = H(i);
            J1 = J(i);
            for(j=1:3)
                n=j;
                G2 = G(j);
                H2 = H(j);
                J2 = J(j);
                eigB(i,j) = (pi^2/sqrt(trho))*( G1^4*Dx/a^4 + G2^4*Dy/b^4 + 2*H1*H2*Dxy/(a^2*b^2) + 4*Dk*(J1*J2 - H1*H2)/(a^2*b^2) )^(1/2);
            end
        end
        eigR=reshape(eigB,[],1);  % reshape 2D array to 1D vector
        eigB=sort(eigR);
        eigB=eigB(1:6);

    case 'ssss'
        for(i=1:3)
            m=i;
            for(j=1:3)
                n=j;
                eigB(i,j) = ( (pi^2/(R^2*b^2))*(1/sqrt(trho)) )*( D(1,1)*(m^4) + 2.*( D(1,2) + 2.*D(3,3) )*(m^2)*(n^2)*(R^2) + D(2,2)*(n^4)*(R^4) )^(1/2);
            end
        end
        eigR=reshape(eigB,[],1);  % reshape 2D array to 1D vector
        eigB=sort(eigR);
        eigB=eigB(1:6);
    case 'cscs'
        for(i=1:3);
            m=i;
            for(j=1:3)
                n=j;
                if( m==1 && n== 1)
                    alpha1 = 4.730;
                    alpha3 = n*pi;
                    alpha2 = 12.30*n^2*pi^2;
                else
                    alpha1 = (m + 0.5)*pi;
                    alpha3 = n*pi;
                    alpha2 = n^2*pi^2*alpha1*(alpha1-2);
                end
                eigB(i,j) = (1/(a^2*sqrt(trho)))*( D(1,1)*(alpha1^4) + 2.*( D(1,2) + 2.*D(3,3) )*(R^2)*alpha2 + D(2,2)*(R^4)*alpha3^4 )^(1/2);
            end
        end
        eigR=reshape(eigB,[],1);  % reshape 2D array to 1D vector
        eigB=sort(eigR);
        eigB=eigB(1:6);

    case 'cccc'
        for(i=1:3);
            m=i;
            for(j=1:3)
                n=j;
                if( m==1 && n== 1)
                    alpha1 = 4.730;
                    alpha3 = 4.730;
                    alpha2 = 151.3;
                else
                    alpha1 = (m + 0.5)*pi;
                    alpha3 = (n + 0.5)*pi;
                    alpha2 = (alpha1*alpha3*(alpha1-2))*(alpha3-2);
                end
                eigB(i,j) = (1/(a^2*sqrt(trho)))*( D(1,1)*(alpha1^4) + 2.*( D(1,2) + 2.*D(3,3) )*(R^2)*alpha2 + D(2,2)*(R^4)*alpha3^4 )^(1/2);
            end
        end
        eigR=reshape(eigB,[],1);  % reshape 2D array to 1D vector
        eigB=sort(eigR);
        eigB=eigB(1:6);

    otherwise
        eigB = [0; 0; 0; 0; 0; 0];
        warning('Unexpected bcType')
end
end
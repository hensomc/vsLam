function [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF2(ipoltype,Nmax,xn,yn)
%% Purpose:  Generate F and its derivative, numerically, vectorized operation

if(ipoltype == 1 || ipoltype == 3)
    [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF_power(Nmax,xn,yn);
end

if(ipoltype == 2)

    if(Nmax >=3 && Nmax <= 8)
%        [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF_legendre3(Nmax,xn,yn);

        % Allocate arrays to save time
        d1=(Nmax+1)^2;
        d2=length(xn);
        W=zeros([d2,d1]);dWx=zeros([d2,d1]);dWy=zeros([d2,d1]);dWxy=zeros([d2,d1]);dWxx=zeros([d2,d1]);dWyy=zeros([d2,d1]);
        
        % Loop over xi,yi
        for i=1:length(xn)
            for j=1:length(yn)
                xi=xn(i); yi=yn(j);
                switch Nmax
                    case 3
                        [W(i,j),dWx(i,j),dWy(i,j),dWxy(i,j),dWxx(i,j),dWyy(i,j)]=FdF_legendre3(Nmax,xi,yi);
                    case 4
                        [W(i,j),dWx(i,j),dWy(i,j),dWxy(i,j),dWxx(i,j),dWyy(i,j)]=FdF_legendre4(Nmax,xi,yi);
                    case 5
                        [W(i,j),dWx(i,j),dWy(i,j),dWxy(i,j),dWxx(i,j),dWyy(i,j)]=FdF_legendre5(Nmax,xi,yi);
                    case 6
                        [W(i,j),dWx(i,j),dWy(i,j),dWxy(i,j),dWxx(i,j),dWyy(i,j)]=FdF_legendre6(Nmax,xi,yi);
                    case 7
                        [W(i,j),dWx(i,j),dWy(i,j),dWxy(i,j),dWxx(i,j),dWyy(i,j)]=FdF_legendre7(Nmax,xi,yi);
                    case 8
                        [W(i,j),dWx(i,j),dWy(i,j),dWxy(i,j),dWxx(i,j),dWyy(i,j)]=FdF_legendre8(Nmax,xi,yi);
                end
            end
        end

%         Option 2 - must restructure FdF_legendre3
%         [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF_legendre3(Nmax,xn,yn);
%         W=reshape(W,d2,d1);dWx=reshape(dWx,d2,d1);dWy=reshape(dWy,d2,d1);dWxy=reshape(dWxy,d2,d1);dWxx=reshape(dWxx,d2,d1);dWyy=reshape(dWyy,d2,d1);

    else
        % Kron approach - slowest
        [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF_legendre(Nmax,xn,yn);
    end
    %disp(['FdF: size dWxx = ',int2str(size(dWxx))])
end
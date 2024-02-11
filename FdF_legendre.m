function [F,dFx,dFy,dFxy,dFxx,dFyy] = FdF_legendre(nDeg, x, y)

    for i=1:length(x)
        [cx, cpx, cdpx] = legendre_poly(nDeg,x(i));
        Px(i,:)=cx;
        dPx(i,:)=cpx;
        dPxx(i,:)=cdpx;
    end

    for i=1:length(y)
        [cy, cpy, cdpy] = legendre_poly(nDeg,y(i));
        Py(i,:)=cy;
        dPy(i,:)=cpy;
        dPyy(i,:)=cdpy;
    end
%     disp('FdF_legendre_poly: Size Px='); size(Px)
%     disp('FdF_legendre_poly: Size Py='); size(Py)

    for i=1:length(x)
        F(i,:)    = kron( Px(i,:), Py(i,:));
        dFx(i,:)  = kron(dPx(i,:), Py(i,:));
        dFy(i,:)  = kron( Px(i,:),dPy(i,:));
        dFxy(i,:) = kron(dPx(i,:),dPy(i,:));
        dFxx(i,:) = kron(dPxx(i,:),Py(i,:));
        dFyy(i,:) = kron( Px(i,:),dPyy(i,:));
        
%         F(i,:)    = kron_fast( Px(i,:), Py(i,:));
%         dFx(i,:)  = kron_fast(dPx(i,:), Py(i,:));
%         dFy(i,:)  = kron_fast( Px(i,:),dPy(i,:));
%         dFxy(i,:) = kron_fast(dPx(i,:),dPy(i,:));
%         dFxx(i,:) = kron_fast(dPxx(i,:),Py(i,:));
%         dFyy(i,:) = kron_fast( Px(i,:),dPyy(i,:));
    end
    
% Attempt to vectorize above step    
%     F    = kron( Px, Py);
%     dFx  = kron(dPx, Py);
%     dFy  = kron( Px,dPy);
%     dFxy = kron(dPx,dPy);
%     dFxx = kron(dPxx,Py);
%     dFyy = kron( Px,dPyy);
    
    % Tensor products
%     F    = kron(Px,Py);
%     dFx  = kron(dPx,Py);
%     dFy  = kron(dPy,Px);
%     dFxy = kron(dPx,dPy);
%     dFxx = kron(dPxx,Py);
%     dFyy = kron(dPyy,Px);
    
%     disp('FdF_legendre_poly: Size F='); size(F)
%     disp('FdF_legendre_poly: Size dFx='); size(dFx)
%     disp('FdF_legendre_poly: Size dFy='); size(dFy)
%     disp('FdF_legendre_poly: Size dFxy='); size(dFxy)
%     disp('FdF_legendre_poly: Size dFxx='); size(dFxx)
%     disp('FdF_legendre_poly: Size dFyy='); size(dFyy)
end
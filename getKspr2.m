function [Kspr]=getKspr2(ipoltype,pDeg,bcType,bcSpring,edgeBC,edgeGID,isoGridID)
Kspr=0;
Kmat=0;
Kt=0;

zvec=zeros(1,(pDeg+1)^2);

for i=1:length(edgeBC)
    % Edge grid IDs and coords
    BC=edgeGID(i,:);
    BCpsieta=[isoGridID(BC,2) isoGridID(BC,3)];
    
    if( bcType(i) ~= 'f' )  % skip if edge is free
        K1=0; K2=0; K3=0;
        for ibc=1:length(BC)
            % Eval FdF using alt method
            psi = BCpsieta(ibc,1);
            eta = BCpsieta(ibc,2);
            
            [W,dWx,dWy,dWxy,dWxx,dWyy]=FdF(ipoltype,pDeg,psi,eta);
            U=W;dUx=dWx;dUy=dWy;dUxy=dWxy;dUxx=dWxx;dUyy=dWyy;
            V=W;dVx=dWx;dVy=dWy;dVxy=dWxy;dVxx=dWxx;dVyy=dWyy;
            P=W;dPx=dWx;dPy=dWy;dPxy=dWxy;dPxx=dWxx;dPyy=dWyy;
            R=W;dRx=dWx;dRy=dWy;dRxy=dWxy;dRxx=dWxx;dRyy=dWyy;
            
            FBC1=double(W);
            FBC2=double(dWx);
            FBC3=double(dWy);
            %FBC = [double(W) double(dWx) double(dWy)];
            
            u=[double(U)    zvec zvec zvec zvec];
            v=[zvec double(V)    zvec zvec zvec];
            w=[zvec zvec    double(W) zvec zvec];
            p=[zvec zvec zvec    double(P) zvec];
            r=[zvec zvec zvec zvec    double(R)];
            
            %Kmat = Kmat + edgeBC(i,:)'*(FBC'*FBC)
            %     Kt = Kt + FBC1'*FBC1;
            K1 = edgeBC(i,1)*(FBC1'*FBC1);
            K2 = edgeBC(i,2)*(FBC2'*FBC2);
            K3 = edgeBC(i,3)*(FBC3'*FBC3);  % NaN for even valued power series
            
            % redefined to use mindlin theory, hardwired to a clamped BC
            K1 = 1.0*(u'*u);
            K2 = 1.0*(v'*v);
            K3 = 1.0*(w'*w);
            K4 = 1.0*(p'*p);
            K5 = 1.0*(r'*r);
            
            %Kt = Kt + K1 + K2 + K3;
            Kt = Kt + K1 + K2 + K3 + K4 + K5;
        end
    end
end
Kspr = bcSpring*Kt;
disp('GenBC: Size of Kspr = '); size(Kspr)

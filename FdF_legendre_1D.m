function [F,dFx,dFxx] = FdF_legendre_1D(nDeg, x)

    %disp('FdF_legendre_1D called');
    
    for i=1:length(x)
        [cx, cpx, cdpx] = legendre_poly(nDeg,x(i));
        Px(i,:)=cx;
        dPx(i,:)=cpx;
        dPxx(i,:)=cdpx;
    end


%     disp('FdF_legendre_poly: Size Px='); size(Px)

    % Form matrices
    F    = Px;
    dFx  = dPx;
    dFxx = dPxx;
end
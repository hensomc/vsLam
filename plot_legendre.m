function plot_legendre(nDeg)
%   getting started
      %clear all, close all,  nfig = 0;
      close all,  nfig = 0;
%
%   set color and marker code for creating plots
      scm = ['r-';      % red solid
             'g:';      % green dotted
             'b-';      % blue solid
             'm:';      % magenta dotted
             'c-';      % cyan solid
             'y:';      % yellow dotted
             'r:';      % red dotted
             'g-';      % greeb solid
             'b:';      % blue dotted
             'y-'];     % yellow solid
%
%   set up independent variable
      Nx = 51;   x = linspace(-1,1,Nx);  y = x; 
      
      % to plot Matlab's legendre
      %for n = 1:nDeg
      %  AP = legendre(n-1,x);
      %  P(:,n) = AP(1,:)';
      %end      
      
      % to plot legendre_poly
      for i=1:length(x)
        [cx, cpx, cdpx] = legendre_poly(nDeg,x(i));
        P(i,:)=cx;
        dPx(i,:)=cpx;    %1st deriv
        dPxx(i,:)=cdpx;  %2nd deriv
      end
      
%   plot curves
      nfig = nfig+1;  figure(nfig)
      for n = 1:nDeg
        plot(x,P(:,n),scm(n,:),'LineWidth',3), grid on, hold on
        txt(n) = {['P',num2str(n-1),'(x)']};  
      end 
      title('\bfFirst 10 Legendre Polynomials','fontsize',14)                               
      xlabel('\bfx Value','fontsize',12),ylabel('\bfP_n(x)','fontsize',12)
      legend(txt,'location','best')
      %legend(txt,'location','northeastoutside')
      
end
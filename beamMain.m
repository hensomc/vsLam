%%  beamMain
% clear workspace, close all plots, clear command window
clear
close all
clc

% Solver inputs
ipol=1;           % 1=power, 2=legendre
pDegI=5;          % initial polynomial degree
pDegF=5;          % final polynomial degree
intMethod = 2;    % 1=simpson integration, 2 = guass quadrature
bcMethod  = 1;    % 1=Explicit Null Space, 2 = springs

plot_flag=false;

% Beam geometry
L=100;    % beam length
b=2;     % beam width
h=2;   % beam height
rho=0.1; % density
E=10E6;   % modulus

% Springs for BCs
k_spring=1.0e8;

i=0;
for pDeg=pDegI:pDegF
    % Solve eigenvalue problem
    [EgN,Ks,Ms,PHI,EG]=beamSolve(ipol,pDeg,intMethod,bcMethod,plot_flag,L,b,h,E,rho,k_spring);
    i=i+1;
    % Save data for graphing
    xval(i,:)=pDeg;
    if( length(EgN) < 5 )
      yval(i,:)=EgN;
    else
      yval(i,1:5)=EgN(1:5);
    end
end

% Convergence plot
if( pDegF ~= pDegI)
    polyType=['Power Series', 'Legendre'];
    figure,plot(xval,yval);
    legend('Mode1','Mode2','Mode3','Mode4','Mode5');
    title('Vibration Frequency Convergence');
    %title('Vibration Frequency Convergence:' + polyType(ipol)); %???
    %title_str='Vibration Frequency Convergence' + polyType(ipol);
    %title(title_str);
    xlabel('Polynomial Degree'),ylabel('Frequency(Hz)');
end


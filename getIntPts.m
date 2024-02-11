function [psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(Mg,intMethod,panXY,plot_flag)

% Define a fine mesh & integration subdivisions
nISubDiv = 20;                                    % number of integration subdivisions

psiIDomain = linspace(-1,1,nISubDiv);             % integration domain x-dir
etaIDomain = linspace(-1,1,nISubDiv);             % integration domain y-dir
dpsi = psiIDomain(2)-psiIDomain(1);                     % dx = distance between 2 integration points
deta = etaIDomain(2)-etaIDomain(1);                     % dy = distance between 2 integration points
[psiIMesh,etaIMesh] = meshgrid(psiIDomain,etaIDomain);

% Nondimensional integration points - default locations, redfined
%                                     if gauss quadrature integration
psiI  = psiIMesh(:) + dpsi/2;        % x-midpoint of integration element
etaI  = etaIMesh(:) + deta/2;        % y-midpoint of integration element

IIDpsi  = find(psiI<1);           % filter out points on boundary psi = 1
psiI   = psiI(IIDpsi);
IIDeta  = find(etaI<1);           % filter out points on boundary eta = 1
etaI   = etaI(IIDeta);

IIDpsieta= intersect(IIDpsi,IIDeta);
psiI   = psiI(IIDpsieta);
etaI   = etaI(IIDpsieta);

gl_wts = 0;


%% Gauss Legendre points and weights
if(intMethod == 2)
    % Uniform mesh
    psiIDomain = linspace(-1,1,nISubDiv);             % integration domain x-dir
    etaIDomain = linspace(-1,1,nISubDiv);             % integration domain y-dir
    [psiIMesh,etaIMesh] = meshgrid(psiIDomain,etaIDomain);
    
    % Gauss integration mesh
    order_1d=[Mg+1,Mg+1];
    order_nd = prod ( order_1d(1:2) );
    [gl_pts,gl_wts]=gl_grid( 2, order_1d, order_nd );
    psiI=gl_pts(1,:);
    etaI=gl_pts(2,:);
end


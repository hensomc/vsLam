%function eigErr = eigErrBuckl(D,a,b,Nx,Ny,EgN,ilam)
function eigErr = eigErrBuckl(panXY,lamDATA,ldsDATA,EgN)

%chk_eigErrBuckl(EgN,ilam)
%    EgN  = Input eigenvalues
%    eigB = Baseline eigenvalues
%    ilam = laminate to be checked (see main driver)

% calculates eigenvalue errors for S-S-S-S orthotropic plate loaded in
% biaxial compression
%

% Inputs:
% D         = D-matrix plate bending rigidities
% m,n       = # buckling waves in x-dir, y-dir
% a,b       = plate dimensions
% Nx,Ny     = in-plane loads
% 
% Outputs:

%recover sol parameters & loads
[ Nx, Ny, Nxy ] = getLdsData( ldsDATA );

abd = get_abd(panXY,0,0,lamDATA);
D=abd(4:6,4:6);

% panel geom
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);

% Get baseline eignvalues
eigB(1,1) = Ncrit_ss(D,1,1,a,b,Nx,Ny)/abs(Nx);
eigB(2,1) = Ncrit_ss(D,2,1,a,b,Nx,Ny)/abs(Nx);
eigB(3,1) = Ncrit_ss(D,3,1,a,b,Nx,Ny)/abs(Nx);
eigB(4,1) = Ncrit_ss(D,2,2,a,b,Nx,Ny)/abs(Nx);
eigB(5,1) = Ncrit_ss(D,1,2,a,b,Nx,Ny)/abs(Nx);
eigB(6,1) = Ncrit_ss(D,1,3,a,b,Nx,Ny)/abs(Nx);
eigB
% Calc error between baseline and input eigenvalues
eigCopy = EgN(1:6)
if(lamDATA.ilam == 1)
   eigErr = abs((eigCopy - eigB)./eigB)*100;
   %eigErr = abs((eigCopy - eigBaseline));
else
   eigErr = [0; 0; 0; 0; 0; 0];
end
%a = 50;
%c = linspace(1,length(eigCopy),length(eigCopy))
%h_eigErr = scatter( eigB, EgN(1:6), a, c, 'filled' );
figure, h_eigErr = scatter( eigB, EgN(1:6));

title( {'\bfCalculated vs Baseline Eigenvalues'; get_prog_vers } )
xlabel('\bfBaseline');
ylabel('\bfCalculated');
grid on;
box on;

% set x&y axes equal
max_eig = max([max(eigB) max(EgN(1:6))]);
min_eig = min([min(eigB) min(EgN(1:6))]);
axis([0, max_eig, 0, max_eig]);

% reference line
refline(1,0);

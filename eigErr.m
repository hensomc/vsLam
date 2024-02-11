function eigError=eigErr(panXY,solDATA,lamDATA,EgN)

%disp('eigError called')

%chk_eigErr(EgN)
%    EgN = Input eigenvalues
%    eigB = Baseline eigenvalues
%    ilam = analysis case to be checked (see main driver)

eigError = [0; 0; 0; 0; 0; 0];

if( lamDATA.ilam > 2 )
    return;
end
%   
ipoltype=solDATA.ipoltype;
pDeg=solDATA.pDeg;
isoltype=solDATA.isoltype;
bcType=solDATA.bcType;

ilam=lamDATA.ilam;
abd=get_abd(panXY, 0, 0, lamDATA);


%% Panel geometry variables
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);

%% get laminate info
[imat, theta0, theta1, thk, mat_names, phi_rot] = readLam(ilam);
[e1,e2,nu12,g12,rho,mtype] = get_mat_props(mat_names);

% build laminate and property strings
nply=length(theta0);
% lamTitle='Laminate: [';
% for i=1:nply
%     if(theta0(i)==theta1(i))
%         thetaStr=num2str(theta0(i));
%     else
%         thetaStr=strcat('<',num2str(theta0(i)),'|',num2str(theta1(i)),'> ');
%     end
%     if(i==1)
%         lamTitle=strcat(lamTitle,thetaStr);
%     else
%         lamTitle=strcat(lamTitle,'/',thetaStr);
%     end
% end
% lamTitle=strcat(lamTitle,']');
% if(e1(1,1)==e2(1,1))
%     propTitle=strcat('E=',num2str(e1(1,1),'%4.1e'), ',  nu=',num2str(nu12(1,1)));
% else
%     propTitle=strcat('E1=',num2str(e1(1,1),'%4.1e'), ',  E2=',num2str(e2(1,1),'%4.1e'), ...
%                      ',  nu12=',num2str(nu12(1,1)),'  G12=',num2str(g12(1,1),'%4.1e'));
% end
% lamPropTitle=strcat(lamTitle, ',   ', propTitle);
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);

%% Get baseline eignvalues
if(ilam == 1)
     eigB = getBaseEig(bcType);
elseif(ilam == 2)
    eigB = getEigOrthoRect(bcType,ilam,a,b,abd);
end

% Copy input eigenvalues to match baseline array size
if(length(EgN) < 6)
    eigB=eigB(1:length(EgN));
    eigCopy(:,1) = EgN(1:length(EgN));
else
    eigCopy(:,1) = EgN(1:6);
end
% if(size(EgN) >= 6)
%   eigCopy = EgN(1:6);
% else
%   eigCopy = EgN;
%   eigB=eigB(1:size(eigCopy), 1);
% end

if(ilam == 1 || ilam == 2)
   eigError = abs((eigCopy - eigB)./eigB)*100;
   %eigErr = abs((eigCopy - eigBaseline));
else
   eigError = [0; 0; 0; 0; 0; 0];
end
%a = 50;
%c = linspace(1,length(eigCopy),length(eigCopy))

%h_eigErr = scatter( eigB, EgN(1:6), a, c, 'filled' );
figure,scatter( eigB, eigCopy);
eigB
box on

%% Titles
title1='\bfRitz Calculated vs Closed Form Eigenvalues';
if(ipoltype == 1 || ipoltype == 3)
    polTitle='Monomial';
elseif(ipoltype == 2)
    polTitle='Legendre';
end
% polTitle = strcat(polTitle, ', N=', num2str(pDeg));
% panTitle =strcat('a=',num2str(a),', b=',num2str(b), ', BC: ',upper(bcType),', ',...
%                   'Polynomial:', polTitle);
% title( {title1;panTitle;lamPropTitle} );
title( {title1;panTitle;lamTitle;polTitle} );
xlabel('\bfReference Frequencies (Hz)');
ylabel('\bfCalculated Frequencies (Hz)');
grid on;

% set x&y axes equal
max_eig = max([max(eigB) max(real(eigCopy))]);
min_eig = min([min(eigB) min(real(eigCopy))]);
axis([0, max_eig, 0, max_eig]);

% reference line
refline(1,0);

    

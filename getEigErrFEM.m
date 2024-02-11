function eigError = getEigErrFEM(panXY, solDATA, lamDATA, nasEigVals ,ritzEigVals ,plot_flag, numEig)

% Initialize 10 eigenvalues to zero
eigError = zeros(10);

isoltype = solDATA.isoltype;

%% Error
%eigError = abs((nasEigVals - EgN(1:10)./EgN(1:10)))*100;

% Table of eigenvalues
disp('===================================')
disp('Eigenvalue comparison: Ritz and FEM');
disp('===================================')
index=min([numEig,length(nasEigVals),length(nasEigVals)]);
eigValRatio = ritzEigVals(1:index)./nasEigVals(1:index);
T1=table( ritzEigVals(1:index), nasEigVals(1:index), eigValRatio(1:index),...
         'VariableNames',{'ritzEigVals','nasEigVals','eigValRatio'});
disp(T1)

%% Scatter graph
figure,scatter( nasEigVals(1:numEig), ritzEigVals(1:numEig));
%figure,scatter( tNasEgN, EgN(1:10));
box on

%% Titles
if(isoltype == 1)
    title1='\bfRitz vs FEM Eigenvalues';
elseif(isoltype == 2)
    title1='\bfRitz vs FEM Buckling Factors';
end
xlabel('\bf\lambda _{FEM}','fontsize',13);
ylabel('\bf\lambda _{RITZ}','fontsize',13);
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);
panLamTitle=strcat(panTitle,', ',lamTitle);
title( {title1;panLamTitle}, 'fontsize',13 );


grid on;

% set x&y axes equal
max_eig = max([max(ritzEigVals(1:numEig)) max(real(nasEigVals(1:numEig) ))]);
min_eig = min([min(ritzEigVals(1:numEig)) min(real(nasEigVals(1:numEig) ))]);
axis([0, max_eig, 0, max_eig]);

% reference line
refline(1,0);

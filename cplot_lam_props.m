function cplot_lam_props( solDATA, panXY, psiBCMesh, etaBCMesh, lamDATA)

% Contour plot of engineering laminate constants
disp('===============CPLOT_LAM_PROPS ===========');

if(solDATA.plot_flag == 0) 
    return
end

% Get titles
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA)
panTitle_lamTitle = strcat(panTitle,',  ',lamTitle)

% Loop over mesh in psi-eta space
% for i=1:length(psiBCMesh)
%     for j=1:length(etaBCMesh)
for i=1:length(etaBCMesh)
    for j=1:length(psiBCMesh)
        abd = get_abd(panXY, psiBCMesh(i,j), etaBCMesh(i,j), lamDATA);
        %abd = get_abd(panXY, -1, 0, lamDATA);
        eij = lam_engr_constants(abd, lamDATA.thk);
        ex(i,j) = eij(1);
        ey(i,j) = eij(2);
        nuxy(i,j) = eij(3);
        gxy(i,j) = eij(4);
    end
end

% Transform mesh
[xBCMesh,yBCMesh] = meshTransform(panXY, psiBCMesh, etaBCMesh);

% Filled contour plot of Ex
figure
pcolor(xBCMesh,yBCMesh,ex);
hold on;
shading interp;
contour(xBCMesh,yBCMesh,ex,'LineStyle','none');
colorbar;

title( {'\bfStiffness Ex';lamTitle;''},'fontsize',12 )
xlabel('\bfx (in.)','fontsize',12)
ylabel('\bfy (in.)','fontsize',12)
hold off;

% Filled contour plot of Ey
figure
pcolor(xBCMesh,yBCMesh,ey);
hold on;
shading interp;
contour(xBCMesh,yBCMesh,ex,'LineStyle','none');
colorbar;

title( {'\bfStiffness Ey';lamTitle;''},'fontsize',12 )
xlabel('\bfx (in.)','fontsize',12)
ylabel('\bfy (in.)','fontsize',12)
hold off;

% Filled contour plot of Gxy
figure
pcolor(xBCMesh,yBCMesh,gxy);
hold on;
shading interp;
contour(xBCMesh,yBCMesh,gxy,'LineStyle','none');
colorbar;

title( {'\bfStiffness Gxy';lamTitle;''},'fontsize',12 )
xlabel('\bfx (in.)','fontsize',12)
ylabel('\bfy (in.)','fontsize',12)
hold off;

end
function getEpsErrFEM(panXY, XYData1, solDATA, lamDATA, nasEps ,ritzEps ,plot_flag)

% Performs comparison of Ritz and FEM strain vectors

isoltype = solDATA.isoltype;

%% Error
%eigError = abs((nasEigVals - EgN(1:10)./EgN(1:10)))*100;

% Table of max/min strains
% format long;
% disp('=============================================')
% disp('Max/Min strain comparison: FEM and Ritz');
% disp('=============================================')
% index=min([numEig,length(nasEigVals),length(nasEigVals)]);
% T1=table( nasEigVals(1:index), ritzEigVals(1:index), ...
%          'VariableNames',{'nasEigVals','ritzEigVals'});

ritzMax = [max(max(ritzEps(:,1))); max(max(ritzEps(:,2))); max(max(ritzEps(:,3))); max(max(ritzEps(:,4))); max(max(ritzEps(:,5)))];
ritzMin = [min(min(ritzEps(:,1))); min(min(ritzEps(:,2))); min(min(ritzEps(:,3))); min(min(ritzEps(:,4))); min(min(ritzEps(:,5)))];
nasMax  = [max(max(nasEps(:,1))); max(max(nasEps(:,2))); max(max(nasEps(:,3))); max(max(nasEps(:,4))); max(max(nasEps(:,5)))];
nasMin  = [min(min(nasEps(:,1))); min(min(nasEps(:,2))); min(min(nasEps(:,3))); min(min(nasEps(:,4))); min(min(nasEps(:,5)))];
rowhdr  = ['Eps_xx'; 'Eps_yy'; 'Eps_xy'; 'Eps_yz'; 'Eps_xz'];

table_title='Max/Min Strain Summary';
disp(table_title);
T1=table (rowhdr, ritzMax, nasMax, ritzMin, nasMin);
format short
disp(T1)

titles={'ritzMax', 'nasMax', 'ritzMin', 'nasMin'};
a=[ritzMax ritzMin nasMax nasMin];
%table2word(titles,a,'vlamWordOutputTables.doc',table_title);

% Section1 cut @ y=b/2
midpt1=0.5*[ panXY(4,1)+panXY(1,1)  panXY(4,2)+panXY(1,2)];
midpt2=0.5*[ panXY(2,1)+panXY(3,1)  panXY(2,2)+panXY(3,2)];
xydata2pt=[1  midpt1;
           2  midpt2];
sectn1=Find_Line_Node([1 2], xydata2pt, XYData1);
xvals=XYData1(sectn1, 2);

% Section2 cut @ x=a/2
midpt1=0.5*[ panXY(1,1)+panXY(2,1)  panXY(1,2)+panXY(2,2)];
midpt2=0.5*[ panXY(3,1)+panXY(4,1)  panXY(3,2)+panXY(4,2)];
xydata2pt=[1  midpt1;
           2  midpt2];
sectn2=Find_Line_Node([1 2], xydata2pt, XYData1);

% Get title strings
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);

% ----------------------
% Graph section1 strains
% ----------------------
% Epsxx
figure
plot( xvals, ritzEps(sectn1,1), 'r:o','linewidth',2), hold on
plot( xvals, nasEps(sectn1,1), 'b-.s','linewidth',2), hold off

title( {'';'\bfTop Fiber Strain \epsilon_{xx} at Midspan as Function of (x)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation x (in.)','fontsize',12)
ylabel('\bf Strain (in/in)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Top Fiber Strain Eps_xx at Midspan as Function of (x)');
T=table (xvals, ritzEps(sectn1,1), nasEps(sectn1,1));
disp(T)


% Epsyy
figure
plot( xvals, ritzEps(sectn1,2), 'r:o','linewidth',2), hold on
plot( xvals, nasEps(sectn1,2), 'b-.s','linewidth',2), hold off

title( {'';'\bfTop Fiber Strain \epsilon_{yy} at Midspan as Function of (x)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation x (in.)','fontsize',12)
ylabel('\bf Strain (in/in)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Top Fiber Strain Eps_yy at Midspan as Function of (x)');
T=table (xvals, ritzEps(sectn1,2), nasEps(sectn1,2));
disp(T)


% Epsxy
figure
plot( xvals, ritzEps(sectn1,3), 'r:o','linewidth',2), hold on
plot( xvals, nasEps(sectn1,3), 'b-.s','linewidth',2), hold off

title( {'';'\bfTop Fiber Strain \epsilon_{xy} at Midspan as Function of (x)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation x (in.)','fontsize',12)
ylabel('\bf Strain (in/in)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Top Fiber Strain Eps_xy at Midspan as Function of (x)');
T=table (xvals, ritzEps(sectn1,3), nasEps(sectn1,3));
disp(T)



% ----------------------
% Graph section2 strains
% ----------------------
% Epsxx
figure
plot( xvals, ritzEps(sectn2,1), 'r:o','linewidth',2), hold on
plot( xvals, nasEps(sectn2,1), 'b-.s','linewidth',2), hold off

title( {'';'\bfTop Fiber Strain \epsilon_{xx} at Midspan as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Strain (in/in)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Top Fiber Strain Eps_xx at Midspan as Function of (y)');
T=table (xvals, ritzEps(sectn2,1), nasEps(sectn2,1));
disp(T)


% Epsyy
figure
plot( xvals, ritzEps(sectn2,2), 'r:o','linewidth',2), hold on
plot( xvals, nasEps(sectn2,2), 'b-.s','linewidth',2), hold off

title( {'';'\bfTop Fiber Strain \epsilon_{yy} at Midspan as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Strain (in/in)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Top Fiber Strain Eps_yy at Midspan as Function of (y)');
T=table (xvals, ritzEps(sectn2,2), nasEps(sectn2,2));
disp(T)

% Epsxy
figure
plot( xvals, ritzEps(sectn2,3), 'r:o','linewidth',2), hold on
plot( xvals, nasEps(sectn2,3), 'b-.s','linewidth',2), hold off

title( {'';'\bfTop Fiber Strain \epsilon_{xy} at Midspan as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Strain (in/in)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Top Fiber Strain Eps_xy at Midspan as Function of (y)');
T=table (xvals, ritzEps(sectn2,3), nasEps(sectn2,3));
disp(T)

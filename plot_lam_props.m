function plot_lam_props( solDATA, panXY, lamDATA)

if(solDATA.plot_flag == 0) 
    return
end

% Get titles
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);
panTitle_lamTitle = strcat(panTitle,',  ',lamTitle);

% rectangular panel case
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);

a=2; b=2;  % use isoparametric panel geom to get ply angles is s-t space
% Calc properties along y=0
x=linspace(-a/2, a/2);
y=zeros(1,length(x));
dx = abs(x(2)-x(1));
ex_eq=0;
for(i=1:length(x))
  abd = get_abd(panXY, x(i), y(i), lamDATA);
  eij = lam_engr_constants(abd, lamDATA.thk);
  ex(i) = eij(1);
  ey(i) = eij(2);
  nuxy(i) = eij(3);
  gxy(i) = eij(4);
  [xc(i),yc(i)] = ISO_st_to_xy(x(i), y(i), panXY );
  ex_eq=ex_eq+ex(i)*dx;
end
ex_eq = ex_eq/a

% Properties as f(x)
h = figure,plot(xc,ex,'r','linewidth',3);hold;
plot(xc,ey,'g','linewidth',3);
plot(xc,gxy,'b','linewidth',3);
%title( {'';'\bfLaminate Engineering Constants as Function of (x)';'';''},'fontsize',14 )
title( {'\bfLaminate Engineering Constants as Function of (x)';panTitle_lamTitle;''},'fontsize',12 )
xlabel('\bfx (in.)','fontsize',12)
ylabel('\bf Stiffness (psi)','fontsize',12)
legend('Ex','Ey','Gxy', 'location','best');

% Properties as f(y)
y=linspace(-b/2, b/2);
x=zeros(1,length(y));
dy = abs(y(2)-y(1));
ey_eq=0;
for(i=1:length(y))
  abd = get_abd(panXY, x(i), y(i), lamDATA);
  eij = lam_engr_constants(abd, lamDATA.thk);
  ex(i) = eij(1);
  ey(i) = eij(2);
  nuxy(i) = eij(3);
  gxy(i) = eij(4);
  [xc(i),yc(i)] = ISO_st_to_xy(x(i), y(i), panXY );
  ey_eq = ey_eq+ey(i)*dy;
end
ey_eq = ey_eq/b;
figure,plot(yc,ex,'r','linewidth',3);hold;
plot(yc,ey,'g','linewidth',3);
plot(yc,gxy,'b','linewidth',3);
%title( {'';'\bfLaminate Engineering Constants as Function of (y)';'';''},'fontsize',14 )
title( {'\bfLaminate Engineering Constants as Function of (y)';panTitle_lamTitle;''},'fontsize',12 )
xlabel('\bfy (in.)','fontsize',12)
ylabel('\bf Stiffness (psi)','fontsize',12)
legend('Ex','Ey','Gxy', 'location','best');

end
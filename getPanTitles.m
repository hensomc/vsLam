function [panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA)

% Recover solution variables from data structure
[ipoltype,pDeg,Mg,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

ilam=lamDATA.ilam;


%% Panel geometry variables
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);

%% get laminate info
%[imat, theta0, theta1, thk, mat_names] = readLam(ilam);
[ilam,imat,theta0,theta1,thk,mat_names,rho_ply,tDeg,thk_coeff]=getLamData(lamDATA);
[e1,e2,nu12,g12,rho,mtype] = get_mat_props(mat_names);

% build laminate and property strings
nply=length(theta0);
%lamTitle='Laminate: [';
lamTitle='Layup: [';

% If nply even, create symmetric laminate title
if(mod(nply,2)==0)
    nply_sym=nply/2;
else
    nply_sym=nply;
end

for i=1:nply_sym
    if(theta0(i)==theta1(i))
        thetaStr=num2str(theta0(i));
    else
        thetaStr=strcat('<',num2str(theta0(i)),'|',num2str(theta1(i)),'> ');
    end
    if(i==1)
        lamTitle=strcat(lamTitle,thetaStr);
    else
        lamTitle=strcat(lamTitle,'/',thetaStr);
    end
end

if(nply_sym == nply)
lamTitle=strcat(lamTitle,']');
else
    lamTitle=strcat(lamTitle,']_{s}');
end

if(e1(1,1)==e2(1,1))
    propTitle=strcat('E=',num2str(e1(1,1),'%4.1e'), ',  nu=',num2str(nu12(1,1)));
else
    propTitle=strcat('E1=',num2str(e1(1,1),'%4.1e'), ',  E2=',num2str(e2(1,1),'%4.1e'), ...
                     ',  nu12=',num2str(nu12(1,1)),'  G12=',num2str(g12(1,1),'%4.1e'));
end
lamPropTitle=strcat(lamTitle, ',   ', propTitle);


%% Titles
if(ipoltype == 1 || ipoltype == 3)
    polTitle='Monomial';
elseif(ipoltype == 2)
    %polTitle='Legendre';
    polTitle='';
end
polTitle = strcat('Ritz:', polTitle, ' Np=', num2str(pDeg), ', Mg=', num2str(Mg));
%panTitle =strcat('a=',num2str(a),', b=',num2str(b), ', BC: ',upper(bcType),', ',polTitle);
panTitle =strcat('BC: ',upper(bcType),',  ',polTitle);


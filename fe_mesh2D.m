function [nasEigVals,nasDisp,nasEps,XYData1]=fe_mesh2D(panXY,solDATA,lamDATA,ldsDATA,nnx,nny,bcType)

disp('Build and run FEM using Nastran...');
%addpath('C:\Users\hensomc\Documents\MATLAB\M_Mesh2D_B');

%% Generate 2D finite element model mesh
keypts = linspace(1,length(panXY),length(panXY));

XYData = [keypts' panXY];

Region = [1 keypts nnx nny];

isoltype = solDATA.isoltype;
% XYData =[     1    -5     0
%      2    20     0
%      3    20    25
%      4    -5    25
%      5     5    10
%      6    15    10
%      7    15    20
%      8     5    20
%      9     5     0
%     10    15     0
%     11    15    25
%     12     5    25];
% Region=[  1 1  9  12  4  3   6
%           2 9  10  6  5   3 3   
%     3    10   2  3  11  3  6
%     4    8  7  11  12  3  2];
% 

%% generate Q4 model: use 4-node quadrilateral elements
ElemType='Q4';
IPLOT=0;
[XYData1,EData1]= M_Mesh2D_B_May7_2015(XYData,Region,ElemType,IPLOT);



%% Display a Table of FE mesh data
%disp('EID   JA   JB   JC    JD  Region')  % header
%disp(EData1)

% Add title to FE Mesh
title('\bfFE mesh: with node and element numbers')

% Label nodes
I=XYData1(:,1);
X=XYData1(I,2);
Y=XYData1(I,3);
text(X,Y,int2str(I),'fontsize',12)

% Label elements
I=EData1(:,1);
JABCD=EData1(I,2:5);
for ie=1:length(I)
    II=I(ie);
    ABCD=JABCD(II,:);
    Xe=sum(X(ABCD))/4;
    Ye=sum(Y(ABCD))/4;
    text(Xe,Ye,int2str(II),'fontsize',12,'color','r')
end

%% Find nodes in regions
NRegion=length(Region(:,1));
for iR=1:NRegion
    nNode=length(XYData1(:,1));
    iE=find(EData1(:,6)==iR);
    ABCD=EData1(iE,2:5);
    AD=ABCD(:);
    In=zeros(nNode,1);
    for i=1:length(AD);
        cc=AD(i)-(1:nNode);
        if min(abs(cc))==0; In(AD(i))=In(AD(i))+1;end
    end
    nodeR{iR}=find(In~=0);
end

% for iR=1:NRegion
%     disp(['Node associated with Region No.',int2str(iR),' are:'])
%     disp(nodeR{iR})
% end

%% Find nodes along region edges - use to apply BCs
edge1_nodes=Find_Line_Node([4 1] ,XYData, XYData1);
edge2_nodes=Find_Line_Node([1 2] ,XYData, XYData1);
edge3_nodes=Find_Line_Node([2 3] ,XYData, XYData1);
edge4_nodes=Find_Line_Node([3 4] ,XYData, XYData1);
edge_BC_nodes = [edge1_nodes; edge2_nodes; edge3_nodes; edge4_nodes];

%% Open file for output
base_name = 'plate_fem';
bdf = [base_name '.bdf'];

fileID=fopen(bdf, 'w');

% SOL card
if( isoltype == 1 ) % eigenvalue
    fprintf(fileID, 'SOL 103\n');
elseif( isoltype == 2 ) % buckling
    fprintf(fileID, 'SOL 105\n');
elseif( isoltype == 3 ) % statics
    fprintf(fileID, 'SOL 101\n');
end

% Slim alter
% fprintf(fileID, '%s\n', 'INCLUDE ''C:\Program Files\TMP\SLIM\SLIM 14.00c\alters\msc\2013.1.0\');
% fprintf(fileID, '%s\n', 'slim105.alter''');

fprintf(fileID, 'CEND\n');

% SUBCASE control
if(isoltype == 1)
    fprintf(fileID, '$-------------------------------------------------------------------------------\n');
    fprintf(fileID, 'SUBCASE 1\n');
    fprintf(fileID, '   METHOD = 1\n');
    fprintf(fileID, '   SPC = 5\n');
    fprintf(fileID, '   WEIGHTCHECK = YES\n');
    
    fprintf(fileID, '   VECTOR(PUNCH,PLOT,SORT1,REAL)=ALL\n');
%    fprintf(fileID, '   SPCFORCES(PLOT,SORT1,REAL)=ALL\n');

elseif(isoltype == 2 || isoltype == 3)
    % Loads/BCs subcase
    fprintf(fileID, '$-------------------------------------------------------------------------------\n');
%     fprintf(fileID, 'MAXMIN(DEF) DISP T1 T2 T3 MAGT R1 R2 R3 MAGR\n');
%     fprintf(fileID, 'MAXMIN(GRID)=ALL\n');
    
    fprintf(fileID, 'SUBCASE 1\n');
    fprintf(fileID, '   LABEL=Load Case 1\n');
    fprintf(fileID, '   SPC = 5\n');
    fprintf(fileID, '   LOAD = 2\n');
    fprintf(fileID, '   WEIGHTCHECK = YES\n');
    fprintf(fileID, '   DISPLACEMENT(PRINT,PUNCH,PLOT,SORT1,REAL)=ALL\n');
        
    %fprintf(fileID, '   DISPLACEMENT(PLOT,SORT1,REAL)=ALL\n');
    fprintf(fileID, '   SPCFORCES(PLOT,SORT1,REAL)=ALL\n');
    %fprintf(fileID, '   STRESS(PLOT,PUNCH,FIBER)=ALL\n');
    %fprintf(fileID, '   STRAIN(PRINT,PLOT,PUNCH,FIBER)=ALL\n');
    %fprintf(fileID, '   GPSTRAIN(PRINT,PLOT,PUNCH,FIBER)=ALL\n');
    %fprintf(fileID, '   STRESS(PRINT,PUNCH,PLOT,FIBER,CORNER)=ALL\n');
    if( isoltype == 3)
        fprintf(fileID, '   STRAIN(PRINT,PLOT,PUNCH,CORNER,FIBER)=ALL\n');
    end
end
if(isoltype == 2)
    % Buckling subcase
    fprintf(fileID, '$-------------------------------------------------------------------------------\n');
    fprintf(fileID, 'SUBCASE 2\n');
    fprintf(fileID, '   METHOD = 1\n');
    fprintf(fileID, '   SPC = 5\n');
    
    fprintf(fileID, '   VECTOR(PLOT,SORT1,REAL)=ALL\n');
    %fprintf(fileID, '   SPCFORCES(PLOT,SORT1,REAL)=ALL\n');

end

fprintf(fileID, '$-------------------------------------------------------------------------------\n');
fprintf(fileID, 'BEGIN BULK\n');

if(isoltype == 1)
    fprintf(fileID, 'PARAM,WTMASS,.002589\n');
elseif(isoltype == 3)
    fprintf(fileID, 'PARAM,NOCOMPS,-1\n');
end

fprintf(fileID, 'PARAM,POST,0\n');   % XDB
%fprintf(fileID, 'PARAM,POST,-1\n');  % OP2

% EIGRL card
if(isoltype == 1 || isoltype == 2)
    fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
    fprintf(fileID, '$EIGRL       SID      V1      V2      ND  MSGLVL  MAXSET  SHFSCL    NORM\n');
    if(isoltype == 1)
        fprintf(fileID, 'EIGRL          1                %8d       0                    MASS\n',solDATA.numMode);
    elseif(isoltype == 2)
        fprintf(fileID, 'EIGRL          1                %8d       0\n',solDATA.numMode);
    end
end
%% Write GRID cards
fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
for i=1:length(XYData1(:,1))
    grid = fprintf(fileID, 'GRID    %8d        %8.3f%8.3f%8.3f\n', XYData1(i,1), X(i), Y(i), 0.0);
    %disp(grid)
end

%% Write CQUAD4 cards
fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
pid=1;
for i=1:length(EData1(:,1))
    if(lamDATA.tDeg > 0 || sum(abs(lamDATA.theta0-lamDATA.theta1))>0)
        pid=EData1(i,1);
    end
    cquad4 = fprintf(fileID, 'CQUAD4  %8d%8d%8d%8d%8d%8d\n',EData1(i,1),pid, EData1(i,2),EData1(i,3),EData1(i,4),EData1(i,5));
    %disp(cquad4)
    
    % Write PSHELL for ea elem if variable thickness
%     if(lamDATA.tDeg > 0)
%         mid=1;
%         JABCD=EData1(i,2:5);
%         Xe=sum(X(JABCD))/4;
%         Ye=sum(Y(JABCD))/4;
%         [psi,eta]=ISO_XY_to_st_Num(Xe,Ye,panXY);
%         t=thk_legendre(lamDATA.tDeg,psi,eta)*lamDATA.thk_coeff';
%         pshell = fprintf(fileID, 'PSHELL  %8d%8d%8.3f%8d\n',pid,mid,t,mid);
%     end
end

%% Write PROPERTY cards
%  Variable thickness
fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
if(lamDATA.tDeg > 0)
    for i=1:length(EData1(:,1))
        pid=EData1(i,1);
        mid=1;
        JABCD=EData1(i,2:5);
        Xe=sum(X(JABCD))/4;
        Ye=sum(Y(JABCD))/4;
        [psi,eta]=ISO_XY_to_st_Num(Xe,Ye,panXY);
        t=thk_legendre(lamDATA.tDeg,psi,eta)*lamDATA.thk_coeff';
        pshell = fprintf(fileID, 'PSHELL  %8d%8d%8.3f%8d\n',pid,mid,t,mid);
    end
end
% Constant thickness
if( lamDATA.tDeg == 0 )
    
    % Recover laminate stack variables
    imat=lamDATA.imat;
    theta0=lamDATA.theta0;
    theta1=lamDATA.theta1;
    thk=lamDATA.thk;
    mat_names=lamDATA.mat_names;
    phi_rot=lamDATA.phi_rot;
    
    nply = length(theta0);
    sout='      NO';
    
    nsm=0.;  

    % Recover panel geometry
%     a=panXY(2,1) - panXY(1,1);
%     b=panXY(3,2) - panXY(2,2);
      a=2; b=2;  % Use isoparametric panel size to calc fiber angles
    
    % Variable fiber path
    if( sum(abs(lamDATA.theta0-lamDATA.theta1)) > 0 )
        for i=1:length(EData1(:,1))
            pid=EData1(i,1);
            JABCD=EData1(i,2:5);
            Xe=sum(X(JABCD))/4;
            Ye=sum(Y(JABCD))/4;
            [psi,eta]=ISO_XY_to_st_Num(Xe,Ye,panXY);

            %fprintf(fileID, 'PCOMP   %8d\n',pid);
            fprintf(fileID, 'PCOMP   %8d        %8.1f\n',pid, nsm);
            for(k=1:nply)
                %theta=get_theta(a,b,theta0(k),theta1(k),phi_rot,psi,eta);
                theta=get_theta(2,2,theta0(k),theta1(k),phi_rot,psi,eta,k);
                fprintf(fileID, '        %8d%8.4f%8.1f%s\n',imat(k),thk(k),theta,sout);
                %disp('theta, psi = '); theta, psi
            end
        end

    else      
        % Single PSHELL
        if(strcmp(lamDATA.mtype(1),'iso') && length(lamDATA.mat_names) == 1)
            pid=1;
            mid=1;
            t=sum(lamDATA.thk);
            %pshell = fprintf(fileID, 'PSHELL  %8d%8d%8.3f%8d%8.3f\n',pid,mid,t,mid,12.*t^3);
            pshell = fprintf(fileID, 'PSHELL  %8d%8d%8.3f%8d\n',pid,mid,t,mid);
            %disp(pshell)
            
        else
            % Single PCOMP
            pid=1;
            nply = length(theta0);
            t=sum(thk);
            z0 = -t/2.;
            sout='      NO';
            
            %fprintf(fileID, 'PCOMP   %8d\n',pid);
            fprintf(fileID, 'PCOMP   %8d        %8.1f\n',pid, nsm);

            for k=1:nply
                fprintf(fileID, '        %8d%8.4f%8.1f%s\n',imat(k),thk(k),theta0(k),sout);
            end
        end
    end
end

%% Write MAT cards
fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
mid=1;
for i=1:length(lamDATA.mat_names);
    if(strcmp(lamDATA.mtype(i),'iso'))
        %mat = fprintf(fileID, 'MAT1    %8d%8.0f%8.0f%8.3f%8.3f\n',mid,lamDATA.e1,lamDATA.g12,lamDATA.nu12,lamDATA.rho);
        mat = fprintf(fileID, 'MAT1    %8d%8.0f        %8.3f%8.3f\n',mid,lamDATA.e1(i),lamDATA.nu12(i),lamDATA.rho(i));
    elseif(strcmp(lamDATA.mtype(i),'ortho'))
        %mat = fprintf(fileID, 'MAT8    %8d%8.0f%8.0f%8.3f%8.0f%8.0f%8.0f%8.3f\n',mid,lamDATA.e1,lamDATA.e2,lamDATA.nu12,lamDATA.g12,lamDATA.g12,lamDATA.g12,lamDATA.rho);
        mat = fprintf(fileID, 'MAT8    %8d%8.0f%8.0f%8.3f%8.0f                %8.3f\n',mid,lamDATA.e1(i),lamDATA.e2(i),lamDATA.nu12(i),lamDATA.g12(i),lamDATA.rho(i));
    end
    %disp(mat)
    mid=mid+1;
end


%% Write SPC cards
% Decode BCs
fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
spcID=0;
spcSets=[];
for(i=1:length(bcType))
    switch bcType(i)
        case {'s'}
%             bcStr = '123';
%             if(isoltype==2 && i==3) % hardwired for Nx-only buckling soln
%                 bcStr='23';
            %bcStr = '3';
%             if(isoltype==2 && i==3) % hardwired for Nx-only buckling soln
%                 bcStr='23';

            %if(i==1) % hardwired for Nx-only buckling soln, or pressure loading
            if(i==2) % hardwired for Ny-only (TEMP)
                bcStr='123';
            else
                bcStr = '3';
            end
        case {'c'}
            bcStr = '123456';
            
            if(isoltype==2)
                %if(i==1)  % hardwired for Nx-only buckling solution, all edges have translation except edge1
                if(i==2)  % hardwired for Ny-only buckling solution (TEMP)
                    bcStr='123456';
                else
                    bcStr='3456';
                end
            end
            
%             if(isoltype==2 && i==3) % hardwired for Nx-only buckling soln- allows edge to move in load direction
%                 bcStr='23456';
%             end
%             if(isoltype==2 && i==4) % hardwired for NY-only buckling soln - allows edge to move in load direction
%                 bcStr='13456';
%             end
            
        case {'f'}
            bcStr = '';
        otherwise
            bcStr = '';
    end
    if( strcmp(bcType(i),'f') )
    else
        spcID=spcID+1;
        spcSets=[spcSets spcID];
        %disp('FEM edge_BC_nodes = '); edge_BC_nodes(i,:)
        %fprintf(fileID, 'SPC1    %8d%8s%8d%8d%8d%8d%8d\n',spcID,bcStr,edge_BC_nodes(i,:));
        %fprintf(fileID, 'SPC1    %8d%8s%8dTHRU%8d\n',spcID,bcStr,edge_BC_nodes(i,1),edge_BC_nodes(i,end));

        if(length(edge_BC_nodes(i,:)) > 6)
            fprintf(fileID, 'SPC1    %8d%8s%8d%8d%8d%8d%8d%8d\n',spcID,bcStr,edge_BC_nodes(i,1:6));
            fprintf(fileID, '        %8d%8d%8d%8d%8d%8d%8d%8d\n',edge_BC_nodes(i,7:end));
            fprintf(fileID, '\n');
        else
            fprintf(fileID, 'SPC1    %8d%8s%8d%8d%8d%8d%8d%8d',spcID,bcStr,edge_BC_nodes(i,1:end));
            fprintf(fileID, '\n');
        end
    end
end
str=strtrim(sprintf('%8d',spcSets));
fprintf(fileID, 'SPCADD         5       %s\n',str);
% fprintf(fileID, 'SPC1           1  123456%8d%8d%8d%8d%8d\n',edge1_nodes(:));
% fprintf(fileID, 'SPCADD         5       1       2       3       4\n');

%% Write LOAD/PRESSURE/FORCE cards
[ Nx, Ny, Nxy, ptLds, ptLoc, p0 ] = getLdsData( ldsDATA );

% Vector of panel edge lengths
panXY3D=[panXY(1,:) 0
         panXY(2,:) 0
         panXY(3,:) 0
         panXY(4,:) 0];

edge_len(1)=norm( panXY(2,:)-panXY(1,:) );
edge_len(2)=norm( panXY(3,:)-panXY(2,:) );
edge_len(3)=norm( panXY(4,:)-panXY(3,:) );
edge_len(4)=norm( panXY(1,:)-panXY(4,:) );

% Element edge length
elem_edge_len(1)=edge_len(1)/(nnx-1);
elem_edge_len(2)=edge_len(2)/(nny-1);
elem_edge_len(3)=edge_len(3)/(nnx-1);
elem_edge_len(4)=edge_len(4)/(nny-1);

% Calculate nodal forces
fx_edge1=Nx*elem_edge_len(2);
fy_edge2=Ny*elem_edge_len(1);
if(abs(Nxy) > 0)
    fxy_edge1=Nxy*elem_edge_len(2);
end

% Nodes with applied loads
if(isoltype == 2) %buckling
    if(abs(Nx) > 0.0)
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$FORCE       SID     GID     CID       F      N1      N2      N3\n');
        % Apply load to edge3_nodes, react with edge 1
        for i=1:length(edge3_nodes)
            if( i==1 || i==length(edge3_nodes))
                fval=0.5*fx_edge1;
            else
                fval=fx_edge1;
            end
            fprintf(fileID, 'FORCE          1%8d        %8.3f%8.3f%8.3f%8.3f\n',edge3_nodes(i),fval,1.,0.,0.);
        end
    end
    
    if(abs(Ny) > 0.0)
        % Apply load to edge4_nodes, react with edge 2
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$FORCE       SID     GID     CID       F      N1      N2      N3\n');
        % Apply load to edge4_nodes, react with edge 2
       for i=1:length(edge4_nodes)
            if( i==1 || i==length(edge4_nodes))
                fval=0.5*fy_edge2;
            else
                fval=fy_edge2;
            end
            fprintf(fileID, 'FORCE          1%8d        %8.3f%8.3f%8.3f%8.3f\n',edge4_nodes(i),fval,0.,1.,0.);
        end
    end
    
    if(abs(Nxy) > 0.0)
    end
    
elseif(isoltype == 3) %static 
    % In-plane loads
    if(abs(Nx) > 0.0)
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$FORCE       SID     GID     CID       F      N1      N2      N3\n');
        % Apply load to edge3_nodes, react with edge 1
        for i=1:length(edge3_nodes)
            if( i==1 || i==length(edge3_nodes))
                fval=0.5*fx_edge1;
            else
                fval=fx_edge1;
            end
            fprintf(fileID, 'FORCE          1%8d        %8.3f%8.3f%8.3f%8.3f\n',edge3_nodes(i),fval,1.,0.,0.);
        end
    end

    if(abs(Ny) > 0.0)
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$FORCE       SID     GID     CID       F      N1      N2      N3\n');
        % Apply load to edge4_nodes, react with edge 2
        for i=1:length(edge4_nodes)
            if( i==1 || i==length(edge4_nodes))
                fval=0.5*fy_edge2;
            else
                fval=fy_edge2;
            end
            fprintf(fileID, 'FORCE          1%8d        %8.3f%8.3f%8.3f%8.3f\n',edge4_nodes(i),fval,0.,1.,0.);
        end
    end
    
    if(abs(Nxy) > 0.0)
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$FORCE       SID     GID     CID       F      N1      N2      N3\n');
        % Apply load to edge3_nodes, react with edge 1
        for i=1:length(edge3_nodes)
            if( i==1 || i==length(edge3_nodes))
                fval=0.5*fxy_edge1;
            else
                fval=fxy_edge1;
            end
            fprintf(fileID, 'FORCE          1%8d        %8.3f%8.3f%8.3f%8.3f\n',edge3_nodes(i),fval,0.,1.,0.);
        end
    end
    
    % Point loads
    if( max(any(abs(ptLds))) > 0 )
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$FORCE       SID     GID     CID       F      N1      N2      N3\n');
        num_ptLds = size(ptLds,1);
        for i=1:num_ptLds
            % search for nodes with pt loads
            %I=XYData1(:,1); X=XYData1(I,2); Y=XYData1(I,3);
            %idx = find(XYData1(I,2) == ptLoc(i,1) && XYData1(I,3) == ptLoc(i,2))
            idx1 = find(XYData1(:,2) == ptLoc(i,1));
            idx2 = find(XYData1(idx1,3) == ptLoc(i,2));
            node_id = XYData1(idx1(idx2),1);
            fprintf(fileID, 'FORCE          1%8d        %8.3f%8.3f%8.3f%8.3f\n',node_id,ptLds(i,3),0.,0.,1.);
        end
    end
    
    % Distributed pressure
    if( max(any(abs(p0))) > 0 )
        fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
        fprintf(fileID, '$PLOAD2      SID       P    EID1    EID2    THRU    EID3\n');
        fprintf(fileID, 'PLOAD2         1%8.3f%8d    THRU%8d\n',p0,I(1),I(end));
    end

    
end

% LOAD CARD
fprintf(fileID, '$---1---|---2---|---3---|---4---|---5---|---6---|---7---|---8---|---9---|---0---|\n');
fprintf(fileID, 'LOAD           2      1.      1.       1\n');

%% ENDDATA
fprintf(fileID, 'ENDDATA\n');
fclose(fileID);

%% Remove previous solutions
f06 = [base_name '.f06'];
f04 = [base_name '.f04'];
op2 = [base_name '.op2'];
xdb = [base_name '.xdb'];
log = [base_name '.log'];
slim= [base_name '.slim'];

warning('off');
delete( f06, f04, op2, xdb, log, slim );
warning('on');

%% Run Nastran
nastran_path='C:\MSC.Software\MSC_Nastran_Student_Edition\2014\Nastran\bin\nastran.exe';
if(exist(nastran_path, 'file') == 0)
    nastran_path='C:\MSC.Software\MSC_Nastran_Student_Edition\20141\Nastran\bin\nastran.exe';
end
cmd = strcat(nastran_path, ' plate_fem.bdf old=no scr=yes 1> nastran_vlam.log');
%nastran_status = system('C:\MSC.Software\MSC_Nastran_Student_Edition\2014\Nastran\bin\nastran.exe plate_fem.bdf old=no scr=yes 1> nastran_vlam.log')
nastran_status = system(cmd);
disp('Nastran job completed, status = ');

%% Recover Nastran Results
[nasEigVals,nasDisp,nasEps]=readNasF06(solDATA, 'plate_fem.f06');

%% Create SLIM DB
slimdb_status = system('"C:\Program Files (x86)\TMP\SLIM\SLIM 14.00c\com\dbgenx" -msc2013.1.0 -i plate_fem.xdb -o plate_fem.slim 1> slim_vlam.log')
disp('Slim DB created');


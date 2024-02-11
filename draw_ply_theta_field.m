function draw_ply_theta_field(panXY,lamDATA,system,numFibPts,refPath,panFibAxis)

disp('Calculating fiber theta field...');

% rectangular panel case
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);
delx=a/10;
dely=b/10;

% Recover laminate data
theta0=lamDATA.theta0;
theta1=lamDATA.theta1;
nply=length(theta0);

for i=1:nply 
    phi_rot=lamDATA.phi_rot;
    k=1; % hardwired to 1st ply
    
    %color='brkgbrkg';
    color='brkggkrbbrkggkrbbrkggkrbbrkggkrb';
    %nply = length(theta0);
    
    [x,y] = meshgrid(panXY(1,1):delx:panXY(2,1), panXY(2,2):dely:panXY(3,2));
    theta = get_theta(a,b,theta0(i),theta1(i),phi_rot,x,y,k)*pi/180;
    u = cos(theta);
    v = sin(theta);
    
    subplot(1,2,1);
    quiver(x,y,u,v);
    pause;
end

%     for(i=1:nply)
%         % if ply is even, change phi_rot to -phi_rot
%         if(mod(i,2) == 0)
%             phi_rot = -phi_rot;
%         end
%             
%         %if(i==1)
%             %h2(i) = draw_fib_paths(theta0(i),theta1(i),phi_rot,color(i),panXY,system, numFibPts,ref_path);
%         %else
%             %check if same ply already drawn
%             drawn=0;
%             for j=1:i-1
%                 if( theta0(i)==theta0(j) && theta1(i)==theta1(j) )
%                     drawn=1;
%                     h2(i)=h2(j);
%                     %disp('ply drawn');
%                     break;
%                 end
%             end
%             if(drawn ==0)
%                 h2(i) = draw_fib_paths(theta0(i),theta1(i),phi_rot,color(i),panXY,system, numFibPts, refPath, panFibAxis);
%             end
%         %end
%         ply_label{i} = strcat('ply-',int2str(i),'<',num2str(theta0(i)),'|',num2str(theta1(i)),'>');
%     end
    disp('FINISHED FIBER THETA FIELD');
    pause;
    
    % Rescale and title x-y fiber paths
    subplot(1,2,1);
    axis( [min(panXY(:,1)) max(panXY(:,1)) min(panXY(:,2)) max(panXY(:,2))] );
    hold on;
    title( {'\bfFiber Paths: (x,y)',' ',' '},'fontsize',14 )

    % Rescale and title s-t fiber paths
    subplot(1,2,2);
    axis( [-1  1 -1  1] );
    hold on;
    %title( {'\bfFiber Paths: (s,t)',' ',' '},'fontsize',14 )
    title( {'\bfFiber Paths: (\xi,\eta)',' ',' '},'fontsize',14 )

end
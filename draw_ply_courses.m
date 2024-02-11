function draw_ply_courses(panXY,lamDATA,system,numFibPts,refPath,panFibAxis)
%draw_ply_courses(panXY,theta0,theta1)

%    panXY = panel corner coords
%    
%        panel origin (0,0)

disp('Calculating fiber paths...');

% rectangular panel case
a=panXY(2,1) - panXY(1,1);
b=panXY(3,2) - panXY(2,2);

% Recover laminate data
theta0=lamDATA.theta0;
theta1=lamDATA.theta1;
phi_rot=lamDATA.phi_rot;

%color='brkgbrkg';
color='brkggkrbbrkggkrbbrkggkrbbrkggkrb';
nply = length(theta0);

if(system == 0)
    % pyhysical panel
    title( {'\bfReference Fiber Paths - Physical System'; get_prog_vers } )
    hold off

elseif(system == 1)
    % equivalent rectangular panel - local panel coordinates
    for(i=1:nply)
        h1(i) = draw_fib_paths(theta0(i),theta1(i),color(i),phi_rot,panXY,system, numFibPts,refPath, panFibAxis);
        ply_label{i} = strcat('ply-',int2str(i),'<',num2str(theta0(i)),'|',num2str(theta1(i)),'>');
    end
    title( {'\bfRef. Fiber Paths - Equiv. Rect. Panel - Local Panel System'; get_prog_vers } )
    legend(h1,ply_label);
    
elseif(system == 2)
    for(i=1:nply)
        % if ply is even, change phi_rot to -phi_rot
        if(mod(i,2) == 0)
            phi_rot = -phi_rot;
        end
            
        %if(i==1)
            %h2(i) = draw_fib_paths(theta0(i),theta1(i),phi_rot,color(i),panXY,system, numFibPts,ref_path);
        %else
            %check if same ply already drawn
            drawn=0;
            for j=1:i-1
                if( theta0(i)==theta0(j) && theta1(i)==theta1(j) )
                    drawn=1;
                    h2(i)=h2(j);
                    %disp('ply drawn');
                    break;
                end
            end
            if(drawn ==0)
                h2(i) = draw_fib_paths(theta0(i),theta1(i),phi_rot,color(i),panXY,system, numFibPts, refPath, panFibAxis);
            end
        %end
        ply_label{i} = strcat('ply-',int2str(i),'<',num2str(theta0(i)),'|',num2str(theta1(i)),'>');
    end
    disp('FINISHED FIBER PATHS');
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
    
    % Legend
    %legend(h2,ply_label,'position','eastoutside');
    %use this: legend(h2,ply_label);

end

end
function h = draw_fib_paths(T0,T1,phi,color,panXY,system,numFibPts,refPath, panFibAxis)
%draw_fiber_paths for a single ply
%    T0    = array of ply angles at panel origin
%    T1    = array of ply angles at panel edges (x=-a/2, a/2)
%    phi   = rigid body rotation of fiber path
%    panXY = panel corner points
%    system = 0,1,2 (rectangular, normalized)
%    refPath = flag to draw only a reference path

% rectangular panel case
%    a = panel length along x-axis
%    b = panel width along y-axis
%    i = y-shift

% Panel length & width
a=abs( ( max( panXY(:,1) ) -min( panXY(:,1) ) ) );
b=abs( ( max( panXY(:,2) ) -min( panXY(:,2) ) ) );

% Panel center
xc=0;yc=0;
for i=1:4
    xc=xc+panXY(i,1);
    yc=yc+panXY(i,2);
end
xc=xc/4;
yc=yc/4;

delx=a/10; %x-shift distance of next fiber path
dely=b/10; %y-shift distance of next fiber path

% Angle to panel point #3 = phi3
phi3 = atand(panXY(3,2)/panXY(3,1));

% ----------------
% theta = theta(y)
% ----------------
% if(T0~=T1 && T0 == 90)
if(panFibAxis == 1)
    %fs=4.0; % fiber spacing
    fs=delx; % fiber spacing
    T0_rot=0;
    %     h=draw_fib_paths2(T0_rot,T1,color,panXY,numFibPts);
    yf=linspace(min(panXY(:,2)),max(panXY(:,2)),numFibPts);
    yf_ref=yf;
    
    for i=0:a
    %for i=-a:a
        % Draw a reference course & reuse it

        % draw in x-y
        if(i==0)
            %xf=linfib3(T0_rot,T1,phi,a,b,yf,i);
            xf=linfib3(T0,T1,phi,a,2*b,yf,i);
        else
            %xf = xf+1;
            xf = xf+delx;
        end
        
        % reposition ref path to panel center    
        if(refPath == 1 && i==0)
            xf_mean=mean(xf);
            xf_ref = xf - xf_mean + xc;
            xf=xf_ref;
        end
        
        subplot(1,2,1)
        hold on;
        
        if( refPath == 1 && i==0 )
            h=plot(xf_ref,yf,color,'LineWidth',3);
        else
            %h=plot(xf,yf,color);
            h=plot(xf_ref-fs*i*cos(phi),yf_ref+fs*i*sin(phi),color);
            h=plot(xf_ref+fs*i*cos(phi),yf_ref-fs*i*sin(phi),color);
        end
        
        if(phi ~= 0)
            [xf,yf]=rot2D( xf, yf, phi, 0., 0. );
        end
        
%         % draw in s-t
%         if(i==0)
%             for k=1:length(yf)
%                 [s(k),t(k)]=ISO_XY_to_st_Num(xf(k),yf(k),panXY);
%             end
%         else
%             s=s+0.2;
%         end
%         
%         subplot(1,2,2)
%         hold on
%         h=plot(s,t,color);
        
        
        %==========
                % draw in s-t
        if(i==0)
            if(T0==T1)  % Straight fibers
                [s1,t1]=ISO_XY_to_st_Num(xf(1),yf(1),panXY);
                [s2,t2]=ISO_XY_to_st_Num(xf(length(xf)),yf(length(yf)),panXY);
                s=linspace(s1,s2,length(xf));
                t=linspace(t1,t2,length(xf));
            else
                for k=1:length(yf)
                    if( refPath == 1)
                        [s_ref(k),t_ref(k)]=ISO_XY_to_st_Num(xf_ref(k),yf(k),panXY);
                        [s(k),t(k)]=ISO_XY_to_st_Num(xf(k),yf(k),panXY);
                    else
                        [s(k),t(k)]=ISO_XY_to_st_Num(xf(k),yf(k),panXY);
                    end
                end
            end
        else
            s=s+0.2;
        end
        
        subplot(1,2,2)
        hold on
        
        if( refPath == 1 && i==0)
            h=plot(s_ref,t_ref,color,'LineWidth',3);
        else
            h=plot(s,t,color);
            h=plot(-s,t,color);
        end
        %==========
        
        
    end
    return
    
elseif(panFibAxis == 0)
    % ----------------
    % theta = theta(x)
    % ----------------
    xf=linspace(min(panXY(:,1)),max(panXY(:,1)),numFibPts);
    xf_ref=xf;
    
    for i=0:b
        %for i=-b:b
        %for i=-b:-b+1
        
        % Draw a reference course & reuse it
        if(i==0)
            
            %Rotate ply course if needed
            if(phi ~= 0)
                T0 = T0 + phi;
                T1 = T1 + phi;
            end
            
            % Get fiber paths
            if(phi <= phi3)
                b=b/cos(phi);
            else
                b=a/sin(phi);
            end
            yf=linfib(T0,T1,phi,a,b,xf,i);
            
            % draw in x-y
            if(T0==T1 && T0==90) % Straight 90-deg fibers
                xf(1:length(xf))=i+a/2;
                yf=linspace(min(panXY(:,2)),max(panXY(:,2)),length(xf));
            end
            
            % Rotate ply course if needed
            %         if(phi ~= 0)
            %             [xf,yf]=rot2D( xf, yf, phi*pi/180., 0., 0. );
            %         end
            
        else
            if(phi ~= 0)
                yf = yf + cos(phi);
                xf = xf - sin(phi);
            else
                %yf = yf+1;
                yf = yf+dely;
            end
        end
        
        %     % Rotate ply course if needed
        %     if(phi ~= 0)
        %         [xf,yf]=rot2D( xf, yf, phi*pi/180., 0., 0. );
        %     end
        
        % reposition ref path to panel center
        if(refPath == 1 && i==0)
            yf_mean=mean(yf);
            yf_ref = yf - yf_mean + yc;
            yf=yf_ref;
        end
        
        subplot(1,2,1)
        hold on
        
        if( refPath == 1 && i==0 )
            h=plot(xf,yf_ref,color,'LineWidth',3);
        else
            %h=plot(xf,yf,color);
            h=plot(xf_ref-i*sin(phi),yf_ref+i*cos(phi),color);
            h=plot(xf_ref+i*sin(phi),yf_ref-i*cos(phi),color);
        end
        
        
        % draw in s-t
        if(i==0)
            if(T0==T1)  % Straight fibers
                [s1,t1]=ISO_XY_to_st_Num(xf(1),yf(1),panXY);
                [s2,t2]=ISO_XY_to_st_Num(xf(length(xf)),yf(length(yf)),panXY);
                s=linspace(s1,s2,length(xf));
                t=linspace(t1,t2,length(xf));
            else
                for k=1:length(xf)
                    if( refPath == 1)
                        [s_ref(k),t_ref(k)]=ISO_XY_to_st_Num(xf(k),yf_ref(k),panXY);
                        [s(k),t(k)]=ISO_XY_to_st_Num(xf(k),yf(k),panXY);
                    else
                        [s(k),t(k)]=ISO_XY_to_st_Num(xf(k),yf(k),panXY);
                    end
                end
            end
        else
            t=t+0.2;
        end
        
        subplot(1,2,2)
        hold on
        
        if( refPath == 1 && i==0)
            h=plot(s_ref,t_ref,color,'LineWidth',3);
        else
            h=plot(s,t,color);
            h=plot(s,-t,color);
        end
        
        %     if(refPath == 1)
        %         return;
        %     end
    end

end


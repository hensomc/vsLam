function theta = get_theta(a,b,T0,T1,phi,x,y,k)
%get_theta(a,b,T0,T1,phi,x,y,k)
%    a = panel length along x-axis
%    b = panel width along y-axis
%    T0 = initial fiber angle
%    T1 = final fiber angle
%    phi= ply rotation
%    k = kth ply index

theta = -9999;

% if ply is even, change phi to -phi
if(mod(k,2) == 0)
    phi = -phi;
end

% ----------------
% theta = theta(y)
% ----------------
% if(T0~=T1 && T0 == 90)
%  
theta=[zeros(size(x))];

for j=1:length(y)
    %     if( (y(j,1) <= b/2) && (0.0 <= y(j,1)) )
    %         theta(j,:) =  2*(abs(T1)-T0)*(y(j,1)/b) + T0;
    %     elseif( (-b/2 <= y(j,1)) && (y(j,1) <= 0.0) )
    %         theta(j,:) = -2*(abs(T1)-T0)*(y(j,1)/b) + T0;
    %     end
    
    %     if( T1 < 0 )
    %         theta(:,j) = -theta;
    %     end
    
    %     if(theta == -9999)
    %         disp('Error calculating theta');
    %     end
    
    theta(j,:)=(T1-T0)*(y(j,1)/b) + T0 + 90;
    
end
return

% end


% ----------------
% theta = theta(x)
% ----------------
    if(T0 == T1)
        theta = T0;
    elseif( x <=(a.*(1.0./2.0) ) && (0.0 <= x) )
            theta = (2.0)*(T1-T0)*(x/a) + T0;
    elseif( (a.*(-1.0./2.0) <= x) && (x <= 0.0) )
            theta = (-2.0)*(T1-T0)*(x/a) + T0;
    end
    
    % rotate ply
    theta = theta + phi;
    
    if(theta == -9999)
        disp('Error calculating theta');
    end
end

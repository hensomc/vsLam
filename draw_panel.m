%function h = draw_panel(panXY,option)
function draw_panel()
%draw_panel(a,b)
% panXY = panel corner point coordinates
%    a  = panel length along x-axis
%    b = panel width along y-axis

%    c = line color
%        panel origin (0,0)

option=1;
panXY = [ 0.0  0.0;  % Square
         10.0  0.0;
         10.0 10.0;
          0.0 10.0];

if(option == 0)
    % draw the panel in physical coordinates
    c1234 = panXY;
    c1    = [panXY(1,1) panXY(1,2)];

elseif(option == 1)
    % rectangular panel: origin (0,0)
    a=panXY(2,1) - panXY(1,1);
    b=panXY(3,2) - panXY(2,2);
    c1 = [a/2, -b/2];
    c2 = [a/2,  b/2];
    c3 = [-a/2, b/2];
    c4 = [-a/2,-b/2];
    c1234 = [c1; c2; c3; c4];
    
elseif(option == 2)
    % normalized coordinates unit square: origin(0,0)
    c1 = [-1, -1];
    c2 = [ 1, -1];
    c3 = [ 1,  1];
    c4 = [-1,  1];
    c1234 = [c1; c2; c3; c4];
end

c12341 = [c1234;c1];
h = figure,plot(c12341(:,1),c12341(:,2));

axis equal; axis tight; hold;


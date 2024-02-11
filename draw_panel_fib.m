function draw_panel_fib(plot_flag, panXY, lamDATA, numFibPts, refPath, panFibAxis)

% panXY = panel corner point coordinates

if(plot_flag > 0)
    %% Fiber Paths
    % Draw a panel in xy & st system
    figure;ISO_PlotQ4_Sub(panXY)
    
    % Draw ply courses in xy and st system
    draw_ply_courses(panXY, lamDATA, 2, numFibPts, refPath, panFibAxis);
    
    % Save to image file
    print('fiber_paths','-dpng');
    
    
    
    
    %% Fiber Theta Field
    figure;ISO_PlotQ4_Sub(panXY)
    
    % Draw ply courses in xy and st system
    draw_ply_theta_field(panXY, lamDATA, numFibPts, refPath, panFibAxis);
    
end
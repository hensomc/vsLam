function gl_grid_display2 (order_1d, points)

%*****************************************************************************80
%
%% GL_GRID_DISPLAY2 displays a 2D Gauss-Legendre grid.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    02 November 2006
%
%  Author:
%
%    John Burkardt
%

%
%  We have to name the axes in order to control the grid.
%
    figure;
    axes_handle = axes;
%
%  Plot the points.
%
    handle = scatter ( points(1,:), points(2,:), 'filled' );
%
%  Force the plotting region to be square, not rectangular.
%
    axis square
%
%  Request grid lines.
%
    grid on
%
%  Specify the location of the grid lines, and suppress labeling.
%
    set ( axes_handle, 'xtick', [ -1, -.75, -.5, -.25, 0, .25, .50, .75, 1] );
    set ( axes_handle, 'xticklabel', [] );
    set ( axes_handle, 'ytick', [ -1, -.75, -.5, -.25, 0, .25, .50, .75, 1] );
    set ( axes_handle, 'yticklabel', [] );
%
%  Make the plotting region slightly bigger than the data.
%
    axis ( [ -1.1, 1.1, -1.1, 1.1 ] )
%
%  Title
%
    s = sprintf ( '%dx%d Gauss-Legendre Grid', order_1d(1), order_1d(2) );
    title ( s );
    
  return
end

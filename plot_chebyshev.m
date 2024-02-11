%function t_polynomial_plot ( index, filename )
function t_polynomial_plot ( nDeg )

%*****************************************************************************80
%
%% T_POLYNOMIAL_PLOT plots Chebyshev polynomials T(n,x).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 March 2012
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer INDEX(*), the orders of 1 or more Chebyshev polynomials
%    to be plotted together.
%
%    Input, string FILENAME, the name into which the graphics information is
%    to be stored.  Note that the PNG format will be used.
%
  a = -1.0;
  b = +1.0;
  m = 501;
  x = linspace ( a, b, m );
  x = x';
  %index_num = length ( index );

  %   set color and marker code for creating plots
      scm = ['r-';      % red solid
             'g:';      % green dotted
             'b-';      % blue solid
             'm:';      % magenta dotted
             'c-';      % cyan solid
             'y:';      % yellow dotted
             'r:';      % red dotted
             'g-';      % greeb solid
             'b:';      % blue dotted
             'y-'];     % yellow solid
%

  clf
  hold on
  %for i = 1 : index_num
  for i = 1 : nDeg
    %n = index(i);
    n=i;
    y = chebyshev_poly ( m, n, x );
    plot ( x, y(:,n+1),scm(n,:), 'LineWidth', 2 );
    txt(n) = {['T',num2str(n-1),'(x)']};  
  end
  grid on
  xlabel ( '<--- X --->' )
  ylabel ( '<--- T(n,x) --->' )
  legend(txt);
  title ( 'Chebyshev polynomials T(n,x)' )
  hold off
  %print ( '-dpng', filename )

  return
end
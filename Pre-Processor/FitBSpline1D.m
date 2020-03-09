function bsp = FitBSpline1D(t, f, n, p)
  
% FitBSpline1D - Fit a B-Spline curve to a set of 1D points.
%
% bsp = FitBSpline1D(t, f, n, p) fit a B-Spline curve to set of given data 
% points (t, f):
%   
%   t   - parametric coordinates of data points                       (in)
%   f   - function values of data points, i.e. f(t)                   (in)
%   n   - number of control points                                    (in)
%   p   - basis degree                                                (in)
%   bsp - struct with B-spline data (n, p, knot, ctrl)                (out)      
%
% See also: EvalBSpline1D, basisfun, findspan

% ==========================================================================
%
%  Copyright (C) 2018  Evandro Parente Junior
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% ==========================================================================
% Created:     16-Jan-2018   Evandro Parente Junior
%
% Modified:
% ==========================================================================

  % Create an open knot vector with uniform increment.
  
  nk = n + p + 1;            % Knot vector size
  ns = n - p;                % Number of knot spans
  dt = (t(end) - t(1))/ns;   % Knot increment
  knot = zeros(1, nk);
  for i = 1:p+1              % Initial knots
    knot(i) = t(1);
  end
  i1 = p + 2;
  i2 = nk - (p + 1);
  for i = i1:i2              % Internal knots
    knot(i) = knot(i-1) + dt;
  end 
  for i = i2+1:nk            % Final knots
    knot(i) = t(end);
  end    
  knot;
  
  % Evaluate the basis functions at the chosen parametric coordinates.

  np = length(t);                        % Number of sampling points 
  s = findspan(n-1, p, t, knot);         % Find the knot span of each t
 
  % Evaluate the basis function matrix [H].

  H = zeros(np, n);
  for i = 1:np
      
    % Get the active B-Spline basis and derivatives.

    B = basisfun(s(i), t(i), p, knot);
      
    % Evaluate the [H] matrix.
    
    ci = s(i) - p;
    for j = 1:p+1
      H(i, ci+j) = B(j);  
    end
  end
  H;
  
  % Evaluate the least square matrix [A].
  
  F = H(:, 2:n-1);
  A = F'*F;
  
  % Evaluate the rhs vector {v}.
  
  v = f' - H(:, 1)*f(1) - H(:, n)*f(np);

  % Evaluate the internal control points solving the normal equations.
  
  intf = A\(F'*v);
  
  % Create the B-Spline.

  bsp.n = n;
  bsp.p = p;
  bsp.knot = knot;
  bsp.ctrl = [f(1) intf' f(np)];
  
  % Evaluate the aproximate values at the given t coordinates.
  
  fa = zeros(1, np);
  for i = 1:np
  
    % Get active control points and B-Spline basis.
    
    ii = s(i) - bsp.p + 1;
    fcp = bsp.ctrl(ii:ii+bsp.p);
    B = basisfun(s(i), t(i), bsp.p, bsp.knot);
    
    % Evaluate the approximate function and the approximation error.
    
    fa(i) = B*fcp';
  end
  fa;
  
  % Evaluate the RMSE and R2 parameters.

  avg = sum(f)/np;
  dif1 = f - fa;
  dif2 = f - avg;
  ssres = sum(dif1.^2);
  sstot = sum(dif2.^2);
  RMSE = ssres/np
  R2 = 1 - ssres/sstot 
end 



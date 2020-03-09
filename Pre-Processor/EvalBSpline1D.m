function fs = EvalBSpline1D(bsp, t)

% EvalBSpline1D - Evaluate a 1D B-Spline function.
%
% fs = EvalBSpline1D(bsp, t) evaluate a 1D B-Spline function at a set of 
% parametric coordinates:
%
%   bsp - struct with B-spline data (n, p, knot, ctrl)                (in)      
%   t   - parametric coordinates                                      (in)
%   fs  - function values, i.e. fs(t)                                 (out)      
%
% See also: FitBSpline1D, basisfun, findspan

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
% Created:     18-Jan-2018   Evandro Parente Junior
%
% Modified:
% ==========================================================================

  % Find the knot span of each t.
 
  s = findspan(bsp.n-1, bsp.p, t, bsp.knot);
  
  % Evaluate each point.
  
  np = length(t);
  fs = zeros(1, np);
  for i = 1:np
  
    % Get active control points and B-Spline basis.
    
    ii = s(i) - bsp.p + 1;
    fcp = bsp.ctrl(ii:ii+bsp.p);
    B = basisfun(s(i), t(i), bsp.p, bsp.knot);
    
    % Evaluate the approximate function and the approximation error.
    
    fs(i) = B*fcp';
  end
end



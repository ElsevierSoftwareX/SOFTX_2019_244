function S = Eshelby_Sphere(nu)

% Eshelby_Sphere - Evaluates the Eshelby tensor for a spherical inclusion.
%
% S = Eshelby_Sphere(nu) evaluates the Eshelby tensor for a spherical 
% inclusion:
%   
%   nu - Poisson's ratio of the matrix phase                         (in)
%   S  - Eshelby tensor in Voigt's notation (6x6 matrix)             (out)      
%
% Ref:  Toshio Mura, "Micromechanics of Defects in Solids", 2nd Edition, 
% Springer Netherlands, pp 588, 1987. ISBN: 978-90-247-3256-2
% DOI: 10.1007/978-94-009-3489-4
%
% See also: Mori_Tanaka

% ==========================================================================
%
%  Copyright (C) 2017  Marcelo Silva Medeiros Junior, Evandro Parente Junior
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
% Created:     29-Aug-2017   Marcelo Silva Medeiros Junior
%
% Modified:    19-Feb-2020   Evandro Parente Junior
%              Code refactoring, comments and formating.
% ==========================================================================

  % Constants depending on Poisson's ratio of the matrix material.
  
  d = (7 - 5*nu)/(15*(1 - nu));
  e = (5*nu - 1)/(15*(1 - nu));
  f = (4 - 5*nu)/(15*(1 - nu));
  
  % Eshelby tensor.
  
  S = [d e e 0 0 0;
       e d e 0 0 0;
       e e d 0 0 0;
       0 0 0 f 0 0;
       0 0 0 0 f 0;
       0 0 0 0 0 f];
end


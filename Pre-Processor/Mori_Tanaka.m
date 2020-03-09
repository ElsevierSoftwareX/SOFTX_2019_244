function C_hom = Mori_Tanaka(C_mat, C_inc, vf, S)

% Mori_Tanaka - Evaluates the homogenized constitutive matrix using the Mori-Tanaka Scheme.
%
% C_hom = Mori_Tanaka(C_mat, C_inc, vf, S) evaluates the elastic  
% constitutive matrix of the composite material using the Mori-Tanaka 
% homogenization scheme:
%   
%   C_mat - constitutive matrix of the matrix phase (6x6)            (in)
%   C_inc - constitutive matrix of the inclusions (6x6)              (in)
%   vf    - volume fraction of the inclusions                        (in)
%   S     - Eshelby tensor in Voigt's notation (6x6 matrix)          (in)      
%   C_hom - constitutive matrix of the composite material (6x6)      (out)     
%
% Ref:  Y. Benveniste, "A New Approach to the Application of Mori-Tanaka's
% Theory in Composite Materials", 1987. DOI: 10.1016/0167-6636(87)90005-6
%
% See also: Eshelby_Sphere, Isotropic_Stiffness, Voigt, Reuss, GSC

% ==========================================================================
%
%  Copyright (C) 2018  Marcelo Silva Medeiros Junior, Evandro Parente Junior
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
% Created:     03-Feb-2018   Marcelo Silva Medeiros Junior
%
% Modified:    06-Mar-2020   Evandro Parente Junior
%              Code refactoring, comments and formating.
% ==========================================================================

  I = eye(6);   % Unit tensor

  % Strain concentration tensor (Eq. 7a).
  
  MT = inv(I + S*inv(C_mat)*(C_inc - C_mat)); 

  % Homogenized moduli (Eq. 14a).
  
  C_hom = C_mat + vf*(C_inc - C_mat)*MT*inv((1 - vf)*I+ vf*MT);
end
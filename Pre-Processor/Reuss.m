function [C_hom] = Reuss(C_mat, C_inc, vf)

% Reuss - Evaluates the homogenized constitutive matrix using the Reuss Scheme.
%
% C_hom = Reuss(C_mat, C_inc, vf) evaluates the elastic  
% constitutive matrix of the composite material using the Voigt 
% homogenization scheme:
%   
%   C_mat - constitutive matrix of the matrix phase (6x6)            (in)
%   C_inc - constitutive matrix of the inclusions (6x6)              (in)
%   vf    - volume fraction of the inclusions                        (in)     
%   C_hom - constitutive matrix of the composite material (6x6)      (out)     
%
% Ref: Jacob Aboudi, Steven M. Arnold, Brett A. Bednarcyk, "Micromechanics 
% of Composite Materials: A Generalized Multiscale Analysis Approach," 
% Butterworth-Heinemann, 2013. ISBN 978-0-12-397035-0.
%
% See also: Isotropic_Stiffness, Voigt, Mori_Tanaka, GSC

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
%              Comments and formating.
% ==========================================================================

  % Apply the Rule of Mixtures to the compliance matrices of the 2 phases.
  
  S_mat = inv(C_mat);
  S_inc = inv(C_inc);
  S_hom = (1 - vf)*S_mat + vf*S_inc;  % Equation (3.66)
  
  % Evaluate the constitutive matrix as the inverse of the compliance matrix.
  
  C_hom = inv(S_hom);
end
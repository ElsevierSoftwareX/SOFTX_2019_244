function [C_hom] = GSC(C_mat, C_inc, vf)

% GSC - Evaluates the homogenized constitutive matrix using the Generalized Self-Consistent Scheme.
%
% C_hom = GSC(C_mat, C_inc, vf) evaluates the elastic constitutive matrix  
% of the composite material using the Generalized Self-Consistent  
% homogenization scheme:
%   
%   C_mat - constitutive matrix of the matrix phase (6x6)            (in)
%   C_inc - constitutive matrix of the inclusions (6x6)              (in)
%   vf    - volume fraction of the inclusions                        (in)      
%   C_hom - constitutive matrix of the composite material (6x6)      (out)     
%
% Ref:  R. Christensen and K. H. Lo, "Solutions for Effective Shear 
% Properties in Three Phase Sphere and Cylinder models", 1979. 
% DOI: 10.1016/0022-5096(79)90032-2
%
% See also: Isotropic_Stiffness, Voigt, Reuss, Mori_Tanaka

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
% Created:     05-Feb-2018   Marcelo Silva Medeiros Junior
%
% Modified:    06-Mar-2020   Evandro Parente Junior
%              Comments and formating.
% ==========================================================================

  % Mechanical properties of matrix and inclusions.
  
  S_mat = inv(C_mat);
  S_inc = inv(C_inc);
  E_m = 1/S_mat(1, 1);
  E_i = 1/S_inc(1, 1);
  nu_m = -E_m*S_mat(1, 2);
  nu_i = -E_i*S_inc(1, 2);
  K_m = E_m/(3 - 6*nu_m);
  K_i = E_i/(3 - 6*nu_i);
  G_m = E_m/(2 + 2*nu_m);
  G_i = E_i/(2 + 2*nu_i);

  % Auxiliary constants - Equation (3.18).

  aux  = G_i/G_m - 1;
  aux1 = G_i/G_m;
  eta1 = (49 - 50*nu_i*nu_m)*aux + 35*aux1*(nu_i - 2*nu_m) + 35*(2*nu_i - nu_m);
  eta2 = 5*nu_i*(aux1 - 8) + 7*(G_i + G_m+4);
  eta3 = aux1*(8 - 10*nu_m) + (7 - 5*nu_m);

  % Coefficients of Equation (3.14): Equations (3.15) - (3.17).

  A = 8*aux*(4-5*nu_m)*eta1*vf^(10/3)-2*(63*aux*eta2+2*eta1*eta3)*vf^(7/3)...
      +252*aux*eta2*vf^(5/3)-50*aux*(7-12*nu_m+8*nu_m^2)*eta2*vf+4*(7-10*nu_m)*eta2*eta3;

  B = -4*aux*(1-5*nu_m)*eta1*vf^(10/3)+4*(63*aux*eta2+2*eta1*eta3)*vf^(7/3)...
      -504*aux*eta2*vf^(5/3)+150*aux*(3-nu_m)*eta2*nu_m*vf+3*(15*nu_m-7)*eta2*eta3;

  C = 4*aux*(5*nu_m-7)*eta1*vf^(10/3)-2*(63*aux*eta2+2*eta1*eta3)*vf^(7/3)...
      +252*aux*eta2*vf^(5/3)+25*aux*(nu_m^2-7)*eta2*vf-(7+5*nu_m)*eta2*eta3;

  % Find the roots of Equation (3.14).

  p = [A B C];
  r = roots(p);
  if (r(1) < 0)
   sol = r(2);
  else 
   sol = r(1);
  end
  
  % Evaluate the mechanical properties of the composite material.
  
  G  = sol*G_m;
  K  = K_m + vf/(1/(K_i - K_m) + (1 - vf)/(K_m + 4*G_m/3));
  E  = (9*K*G)/(3*K + G);
  nu = (3*K - 2*G)/(2*(3*K + G));

  % Evaluate the elastic matrix of the composite material.
  
  C_hom = Isotropic_Stiffness(E, nu);
end
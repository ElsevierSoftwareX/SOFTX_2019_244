
% ==========================================================================
% DemoMicroFEA.m - file to demonstrate the use of homonization schemes 
% available in MicroFEA 1.0.
% ==========================================================================
%
%  Copyright (C) 2020  Marcelo Silva Medeiros Junior, Evandro Parente Junior
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
% Created:     06-Mar-2020   Marcelo Medeiros/Evandro Parente
%
% Modified:    
% ==========================================================================

clear all;
clc;

% Material properties (matrix = Aluminum and inclusions = Al2O3).

E_mat = 70e9;
E_inc = 380e9;
nu_mat = 0.3;
nu_inc = 0.3;

% Elastic matrices.

C_mat = Isotropic_Stiffness(E_mat, nu_mat);
C_inc = Isotropic_Stiffness(E_inc, nu_inc);

% Eshelby tensor.

H = Eshelby_Sphere(nu_mat);

% Volume fractions.

nf = 31;
vf = linspace(0, 1, nf);

% Evaluate the composite properties by homogenization.

E  = zeros(4, nf);
nu = zeros(4, nf);
for j = 1:nf
  C = Voigt(C_mat, C_inc, vf(j));
  S = inv(C);
  E(1, j) = 1/S(1, 1);
  nu(1, j) = -E(1, j)*S(1, 2);
  
  C = Reuss(C_mat, C_inc, vf(j));
  S = inv(C);
  E(2, j) = 1/S(1, 1);
  nu(2, j) = -E(2, j)*S(1, 2);

  C = Mori_Tanaka(C_mat, C_inc, vf(j), H);
  S = inv(C);
  E(3, j) = 1/S(1, 1);
  nu(3, j) = -E(3, j)*S(1, 2);  
  
  C = GSC(C_mat, C_inc, vf(j));
  S = inv(C);
  E(4, j) = 1/S(1, 1);
  nu(4, j) = -E(4, j)*S(1, 2);  
end
E;
nu;

% Evaluate the bulk and shear modulus.

K = E./(3*(1 - 2*nu));
G = E./(2*(1 + nu));

% Plot the composite properties.

set(gcf, 'Name', 'MicroFEA 1.0');
subplot(2, 2, 1);
plot(vf, E, 'LineWidth', 2);
legend('Voigt', 'Reuss', 'Mori-Tanaka', 'Gen. Self-Consistent', 'Location', 'Best');
xlabel('Vf');
ylabel('E (Pa)');
grid on;

subplot(2, 2, 2);
plot(vf, nu, 'LineWidth', 2);
legend('Voigt', 'Reuss', 'Mori-Tanaka', 'Gen. Self-Consistent', 'Location', 'Best');
xlabel('Vf');
ylabel('nu');
grid on;

subplot(2, 2, 3);
plot(vf, K, 'LineWidth', 2);
legend('Voigt', 'Reuss', 'Mori-Tanaka', 'Gen. Self-Consistent', 'Location', 'Best');
xlabel('Vf');
ylabel('K (Pa)');
grid on;

subplot(2, 2, 4);
plot(vf, G, 'LineWidth', 2);
legend('Voigt', 'Reuss', 'Mori-Tanaka', 'Gen. Self-Consistent', 'Location', 'Best');
xlabel('Vf');
ylabel('G (Pa)');
grid on;



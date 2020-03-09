% ==========================================================================
% FitWatanabe2015.m - Fit the experimental volume fraction variation using a 
% B-Spline function.
% ==========================================================================
% Ref: Y. Watanabe, P. D. Sequeira, H. Sato, T. Inamura, H. Hosoda, 
% "Aluminum matrix texture in al{al3ti functionally graded materials 
% analyzed by electron back-scattering diffraction, Japanese Journal of 
% Applied Physics, 2015. DOI:10.7567/jjap.55.01ag03
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

% Experimental data (parametric coordinate x volume fraction).

t  = [0.0000  0.150  0.250  0.350  0.450  0.5500  0.650  0.750  0.850 1];
Vf = [0.0276  0.016  0.016  0.030  0.058  0.0451  0.134  0.155  0.209 0.288];

% Fit a B-Spline to the experimental data.

n  = 7;                          % Number of control points
p  = 3;                          % Basis degree
bsp = FitBSpline1D(t, Vf, n, p)  % B-Spline fit

% Plot the experimental data and the fitted curve.
  
ts = linspace(min(t), max(t), 31);
fs = EvalBSpline1D(bsp, ts);
plot(t, Vf, ts, fs);
xlabel('t');
ylabel('Vf');
legend('Experimental', 'B-Spline', 'Location', 'Best');
title('Volume fraction');
grid on;

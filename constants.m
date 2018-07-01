% Copyright (C) 2017  Max Gattringer

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [resonator] = constants(varargin);
JAN = 0;
nVarargs = length(varargin);
for k = 1:nVarargs
    switch varargin{k}
        case 'janezic' 
            JAN = 1;
    end
end
% RESONATOR DIMENSIONS FOR QUICK MEASUREMENTS

% Please note that for accurate measurements multiple measurements are
% necessary. The results of a single measurement can deviate strongly from
% a "true" value derived from multiple measurements

% Resonator dimensions
resonator.a_u = 38.18075e-3/2;%0.019090375000000

resonator.b = 35e-3; % Flange width

resonator.L_u = 25.009e-3;
resonator.L_l = 25.037e-3;

% Number of modes
resonator.N = 75;

% Calibration % Corrected Calibration!!!!!!
resonator.a_l= 0.019054912915683; % According to RepTest Calibration 
resonator.sigma = 1.005331660728065e+07;



if(JAN)
    L_janezic =  mean([resonator.L_u resonator.L_l]);
    a_janezic = 0.019072613028406; % According to RepTest Calibration
    sigma_janezic = 1.005411529033119e+07;
    
    % Compatability
    resonator.L_u = L_janezic;
    resonator.L_l = L_janezic;
    resonator.sigma = sigma_janezic;
    
    resonator.a_l = a_janezic;
    resonator.a_u = a_janezic;
end

resonator.sigma_f = resonator.sigma; % Old idea to model flange loss separately
end
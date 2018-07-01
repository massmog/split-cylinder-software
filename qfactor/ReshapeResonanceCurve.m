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

function [ f_out, S21_out ] = ReshapeResonanceCurve( f, S21, Delta)
    % Reshape resonance curve around f_r - Using 3dB Method as a first estimate
    [f_r, Q]=ThreeDBmethod(f, S21);

    d_f = f_r*(Delta/Q);

    [min_l, l_index] = min(abs(f-f_r+d_f));
    [min_u, u_index] = min(abs(f-f_r-d_f));
    
    f_out = f(l_index:u_index);
    S21_out = S21(l_index:u_index);
end


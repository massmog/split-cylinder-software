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

function [ f_r, Q, S_max ] = Estin( f_in, S_in )
    % Check input
    if(~isvector(f_in)&&~isvector(S_in))
        error('Input variables must be vectors!');
    end
    % Transform row vectors to column vectors
    if(~iscolumn(f_in))
        f_in = f_in.';
    end
    
    if(~iscolumn(S_in))
        S_in = S_in.';        
    end
        
    % ESTIN Estimator of Q and f_r

    % Reshape resonance curve for correct LS linear fit
    % Frequency spacing factor of Delta = 2.6 (Not optimized!)
    Delta   = 2.9;
    [f,S]   = ReshapeResonanceCurve(f_in, S_in, Delta);
    
    % Find maximum of resonance curve
    [S_max,x_0] = max(abs(S).^2);
    f_r     = f(x_0);
    % Evaluate linearization expression (model y_k = Q.x_k)
    x_k     = abs(f./f_r - f_r./f);
    y_k     = real(sqrt(S_max./(abs(S).^2)-1));
    % Linear Least-Square Estimator
    Q       = (x_k.'*x_k)^-1*x_k.'*y_k;
end


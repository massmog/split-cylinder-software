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

function [ f_r, Q] = ThreeDBmethod( f_in, S_in )
    % Check input
    if(~isvector(f_in) && ~isvector(S_in))
        error('Input variables must be vectors!');
    end
    % Transform row vectors to column vectors
    if(~iscolumn(f_in))
        f_in = f_in.';
    end
    
    if(~iscolumn(S_in))
        S_in = S_in.';        
    end
    % Traditional 3dB method to determine f_r and Q
    % Find resonant frequency
    [S_max,index]   = max(abs(S_in).^2);
    f_r             = f_in(index);

    % Bandwidth algorithm
    f_l             = find(S_max/2 > abs(S_in(1:index)).^2,1,'last');
    f_h             = index+find(S_max/2 < abs(S_in(index:end)).^2,1,'last')-1;

    if(~isempty(f_l)&&~isempty(f_h))
        % Calculate Q
        % This bandwidth estimation algorithm should allow to find a good
        % estimate the Q even for disrupted resonance curves
        BW_est = 2*min([f_in(f_h)-f_r f_r-f_in(f_l)]);
        if(BW_est==0)
            BW_est = f_in(f_h)-f_in(f_l);
        end        
        BW          = BW_est;
        Q           = f_r/BW;
    else
        warning('Estimating Q failed!');
        Q           = 0;
    end

end


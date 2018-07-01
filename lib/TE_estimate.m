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

function [ out ] = TE_estimate( n, m, resonator, f_r, d )
if(n==0)
    temp = zerobess('DJ',n,m+1)';
    pnmd = temp(end);
else
    temp = zerobess('DJ',n,m)';
    pnmd = temp(end);
end

c_0 = 299792458;

L = mean([resonator.L_l,resonator.L_u]);
a = mean([resonator.a_l,resonator.a_u]);

omega = 2*pi*f_r;

res_cond = @(e_r) sqrt(e_r.*(omega./c_0).^2-(pnmd./a).^2).*tan(d./2.*sqrt(e_r.*(omega./c_0).^2-(pnmd./a).^2))-sqrt((omega./c_0).^2-(pnmd./a).^2).*cot(L.*sqrt((omega./c_0).^2-(pnmd./a).^2));
out = fzero(res_cond,10);

end


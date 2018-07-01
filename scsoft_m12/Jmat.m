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

function [ Z, coefficients] = Jmat(resonator, e_r, f, d, varargin)
cond = 1;

a_u = resonator.a_u;
L_u = resonator.L_u;

%% PARAMETER OVERRIDE FUNCTION
nVarargs = length(varargin);
for k = 1:nVarargs
    switch varargin{k} 
        case 'a_u'
            a_u = varargin{k+1};
        case 'L_u'
            L_u = varargin{k+1};
        case 'N_s'
            N_s = varargin{k+1};
        case 'cond_off'
            cond = 0;
    end
end

b = resonator.b;
N = resonator.N;

omega = 2.*pi.*f;
c_0 = 299792458;
e_lab = 1.00055;

k_u = omega./c_0.*sqrt(e_lab);
k_s = omega./c_0.*sqrt(e_r);

N_u = N;
h_n_u = j1(1:N_u)./a_u;
p_n_u = (sqrt(k_u.^2.-h_n_u.^2));

% Compute number of modes from ideal mode ratio
% Optimise number of modes for ideal relative convergence
if(~exist('N_s','var'))
    [d_p,n_modes]=min(abs(imag((sqrt(k_s.^2.-(j1(1:250)./b).^2)))-imag(p_n_u(N_u))));
    if(d_p<0)
        n_modes = n_modes + 1;
    end
    N_s = n_modes;
end

% Set optimum modes
h_m_s = j1(1:N_s)./b;
p_m_s = (sqrt(k_s.^2.-h_m_s.^2)); 

% Constant if row-scaling not used 
if(cond)
    U_n = (p_n_u(N_u))./(cosh(imag(p_n_u).*L_u));
    V_m = (p_m_s(N_s))./(cosh(imag(p_m_s).*d./2));
else
    U_n = 1;
    V_m = 1;
end

R=diag(V_m.*b.^2./2.*besselj(0,h_m_s.*b).^2.*cos(p_m_s.*d./2));

S=diag(U_n.*(p_n_u).*a_u.^2./2.*besselj(0,h_n_u.*a_u).^2.*cos(p_n_u.*L_u));

Q = (a_u.*ones(1,N_s)'*h_n_u)./((h_m_s'*ones(1,N_u)).^2.-(ones(1,N_s)'*h_n_u).^2).*(besselj(1,h_m_s'.*a_u)*(U_n.*besselj(0,h_n_u.*a_u).*sin(p_n_u.*L_u)));
P = (a_u.*h_n_u'*(p_m_s))./((ones(1,N_u)'*h_m_s).^2.-(h_n_u'*ones(1,N_s)).^2).*(besselj(0,h_n_u'.*a_u)*(V_m.*besselj(1,h_m_s.*a_u).*sin(p_m_s.*d./2)));

Z=[Q -R ; S -P];

coefficients.U_n = U_n;
coefficients.h_n_u = h_n_u;
coefficients.p_n_u = p_n_u;

coefficients.V_m = V_m;
coefficients.h_m_s = h_m_s;
coefficients.p_m_s = p_m_s;

coefficients.f = f;
coefficients.e_r = e_r;
coefficients.d = d;

coefficients.N_u = N_u;
coefficients.N_s = N_s;
end


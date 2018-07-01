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

function [ Z, coefficients] = Zmat(resonator, e_r, f, d, varargin)
a_l = resonator.a_l;
a_u = resonator.a_u;
L_l = resonator.L_l;
L_u = resonator.L_u;

%% PARAMETER OVERRIDE FUNCTION
nVarargs = length(varargin);
for k = 1:nVarargs
    switch varargin{k}
        case 'a_l'
            a_l = varargin{k+1};   
        case 'a_u'
            a_u = varargin{k+1};
        case 'L_l'
            L_l = varargin{k+1};
        case 'L_u'
            L_u = varargin{k+1};
        case 'N_s'
            N_s = varargin{k+1};
    end
end

b = resonator.b;
N = resonator.N;

omega = 2.*pi.*f;
c_0 = 299792458;
e_lab = 1.00055;

% Compute number of modes from ideal mode ratio
N_u = N;
N_l = N;

h_n_u = j1(1:N_u)./a_u;
h_n_l = j1(1:N_l)./a_l;

k_u = omega./c_0.*sqrt(e_lab);
k_s = omega./c_0.*sqrt(e_r);

p_n_u = (sqrt(k_u.^2.-h_n_u.^2));
p_n_l = (sqrt(k_u.^2.-h_n_l.^2));

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
U_n = (p_n_u(N_u))./(cosh(imag(p_n_u).*L_u));
V_m = (p_m_s(N_s))./(cosh(imag(p_m_s).*d./2));
W_m = (p_m_s(N_s))./(cosh(imag(p_m_s).*d./2));
L_n = (p_n_l(N_l))./(cosh(imag(p_n_l).*L_l));

M_1 = (a_u.*ones(1,N_s)'*h_n_u)./((h_m_s'*ones(1,N_u)).^2.-(ones(1,N_s)'*h_n_u).^2).*(besselj(1,h_m_s'.*a_u)*(U_n.*besselj(0,h_n_u.*a_u).*sin(p_n_u.*L_u)));
M_2 = diag(V_m.*b.^2./2.*besselj(0,h_m_s.*b).^2.*cos(p_m_s.*d./2));
M_3 = diag(W_m.*b.^2./2.*besselj(0,h_m_s.*b).^2.*sin(p_m_s.*d./2));
M_4 = M_2;
M_5 = -M_3;
M_6 = (a_l.*ones(1,N_s)'*h_n_l)./((h_m_s'*ones(1,N_l)).^2.-(ones(1,N_s)'*h_n_l).^2).*(besselj(1,h_m_s'.*a_l)*(L_n.*besselj(0,h_n_l.*a_l).*sin(p_n_l.*L_l)));
M_7 = diag(U_n.*(p_n_u).*a_u.^2./2.*besselj(0,h_n_u.*a_u).^2.*cos(p_n_u.*L_u));
M_8 = (a_u.*h_n_u'*(p_m_s))./((ones(1,N_u)'*h_m_s).^2.-(h_n_u'*ones(1,N_s)).^2).*(besselj(0,h_n_u'.*a_u)*(V_m.*besselj(1,h_m_s.*a_u).*sin(p_m_s.*d./2)));
M_9 = -(a_u.*h_n_u'*(p_m_s))./((ones(1,N_u)'*h_m_s).^2.-(h_n_u'*ones(1,N_s)).^2).*(besselj(0,h_n_u'.*a_u)*(W_m.*besselj(1,h_m_s.*a_u).*cos(p_m_s.*d./2)));
M_10 = (a_l.*h_n_l'*(p_m_s))./((ones(1,N_l)'*h_m_s).^2.-(h_n_l'*ones(1,N_s)).^2).*(besselj(0,h_n_l'.*a_l)*(V_m.*besselj(1,h_m_s.*a_l).*sin(p_m_s.*d./2)));
M_11 = (a_l.*h_n_l'*(p_m_s))./((ones(1,N_l)'*h_m_s).^2.-(h_n_l'*ones(1,N_s)).^2).*(besselj(0,h_n_l'.*a_l)*(W_m.*besselj(1,h_m_s.*a_l).*cos(p_m_s.*d./2)));
M_12 = diag(L_n.*(p_n_l).*a_l.^2./2.*besselj(0,h_n_l.*a_l).^2.*cos(p_n_l.*L_l));

Z=[M_1,             -M_2,   -M_3,   zeros(N_s,N_l);...
   zeros(N_s,N_u),  -M_4,   -M_5,   M_6;...
   M_7,             -M_8,   -M_9,   zeros(N_u,N_l);...
   zeros(N_l,N_u),  -M_10,  -M_11,  M_12];

coefficients.U_n = U_n;
coefficients.h_n_u = h_n_u;
coefficients.p_n_u = p_n_u;

coefficients.V_m = V_m;
coefficients.W_m = W_m;
coefficients.h_m_s = h_m_s;
coefficients.p_m_s = p_m_s;

coefficients.L_n = L_n;
coefficients.h_n_l = h_n_l;
coefficients.p_n_l = p_n_l;

coefficients.f = f;
coefficients.e_r = e_r;
coefficients.d = d;

coefficients.N_u = N_u;
coefficients.N_s = N_s;
end


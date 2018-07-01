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

function [ Z, coefficients] = Cmat(resonator, f, varargin)
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
            L_u = varargin{k+1};
        case 'L_u'
            L_l = varargin{k+1};
    end
end

N_u = resonator.N;
N_l = resonator.N;
N_s = 0;

omega = 2.*pi.*f;
c_0 = 299792458;
e_lab = 1.00055;

h_n_u = j1(1:N_u)./a_u;
h_n_l = j1(1:N_l)./a_l;

k_u = omega./c_0.*sqrt(e_lab);

p_n_u = (sqrt(k_u.^2.-h_n_u.^2));
p_n_l = (sqrt(k_u.^2.-h_n_l.^2));

% Constant if row-scaling not used 
U_n = (p_n_u(N_u))./(cosh(imag(p_n_u).*L_u));
L_n = (p_n_l(N_l))./(cosh(imag(p_n_l).*L_l));

C_1 = diag(U_n.*a_u.^2./2.*besselj(0,h_n_u.*a_u).^2.*sin(p_n_u.*L_u));
C_2 = -(a_l.*ones(1,N_u)'*h_n_l)./((ones(1,N_u)'*h_n_l).^2.-(h_n_u'*ones(1,N_l)).^2).*(besselj(1,h_n_u'.*a_l)*(L_n.*besselj(0,h_n_l.*a_l).*sin(p_n_l.*L_l)));

C_3 = (a_l*h_n_l.'*p_n_u)./((h_n_l'*ones(1,N_u)).^2.-(ones(1,N_l)'*h_n_u).^2).*((besselj(0,h_n_l'.*a_l))*(U_n.*besselj(1,h_n_u.*a_l).*cos(p_n_u.*L_u)));
C_4 = diag(L_n.*a_l.^2./2.*besselj(0,h_n_l.*a_l).^2.*p_n_l.*cos(p_n_l.*L_l));

Z=[C_1,-C_2;...
   C_3,-C_4];

coefficients.f = f;
coefficients.d = 0;

coefficients.U_n = U_n;
coefficients.h_n_u = h_n_u;
coefficients.p_n_u = p_n_u;

coefficients.L_n = L_n;
coefficients.h_n_l = h_n_l;
coefficients.p_n_l = p_n_l;

coefficients.N_u = N_u;
coefficients.N_l = N_l;
end


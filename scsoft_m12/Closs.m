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

function [ out ] = Closs( resonator, coefficients, null_Z, varargin)
a_l = resonator.a_l;
a_u = resonator.a_u;
L_u = resonator.L_u;
L_l = resonator.L_l;
f_r = coefficients.f;

%% INTEGRATION SETTING
int_atol = 1e-10;
int_rtol = 1e-6;
% Analytical integration switch
ANAL_INT = 1;
%% PARAMETER OVERRIDE FUNCTION
nVarargs = length(varargin);
for k = 1:nVarargs
    switch varargin{k}
        case 'Q' 
            MODE = 1; %% Q MODE
            sigma = resonator.sigma;
            R_s = sqrt(omega*m_0/(2*sigma));
        case 'sigma'
            MODE = 2; %% SIGMA MODE
            Q = varargin{k+1};
        case 'a_l'
            a_l = varargin{k+1};   
        case 'a_u'
            a_u = varargin{k+1};
        case 'L_l'
            L_l = varargin{k+1};
        case 'L_u'
            L_u = varargin{k+1};
        case 'f'
            f_r = varargin{k+1};
    end
end

% Coefficients
d = coefficients.d;

U_n = coefficients.U_n;
h_n_u = coefficients.h_n_u;
p_n_u = coefficients.p_n_u;

L_n = coefficients.L_n;
h_n_l =  coefficients.h_n_l;
p_n_l = coefficients.p_n_l;

% Function tanD - Calculate tangent delta from mode matching problem
b = resonator.b;

N_u = coefficients.N_u;
N_l = coefficients.N_u;

A_n = null_Z(1:N_u);
D_n = null_Z(N_u+1:N_u+N_u);

omega = 2.*pi.*f_r;
e_0 = 8.85418782e-12;
e_lab = 1.00055;
m_0 = 4*pi*1e-7;

% Energy in the Upper Cavity
if(ANAL_INT)
    W_u = sum(e_0*e_lab*pi*a_u^2/4.*(abs(A_n.*U_n').*besselj(0,h_n_u'.*a_u)).^2.*((L_u/2)*ones(N_u,1)-sin(2.*p_n_u.'.*L_u)./(4.*p_n_u.')).*(-(imag(p_n_u')~=0)+(imag(p_n_u')==0)),1);
else
    E_u = @(r,z) sum(A_n.*U_n.'.*besselj(1,h_n_u.'.*r).*sin(p_n_u.'.*(L_u+d/2-z)),1);
    E_u_integrand = @(r1,z1) arrayfun(@(r,z) r.*abs(E_u(r,z)).^2,r1,z1);
    W_u = e_0*e_lab./4.*2.*pi.*integral2(E_u_integrand,0,a_u,d/2,L_u+d/2,'RelTol',int_rtol,'AbsTol',int_atol); 
end
% Energy in the Lower Cavity
if(ANAL_INT)
    W_l = sum(e_0*e_lab*pi*a_l^2/4.*(abs(D_n.*L_n').*besselj(0,h_n_l'.*a_l)).^2.*((L_l/2)*ones(N_u,1)-sin(2.*p_n_l.'.*L_l)./(4.*p_n_l.')).*(-(imag(p_n_l')~=0)+(imag(p_n_l')==0)),1);
else
    E_l = @(r,z) sum(D_n.*L_n.'.*besselj(1,h_n_l.'.*r).*sin(p_n_l.'.*(L_l+d/2+z)),1);
    E_l_integrand = @(r1,z1) arrayfun(@(r,z) r.*abs(E_l(r,z)).^2,r1,z1);
    W_l = e_0*e_lab./4.*2.*pi.*integral2(E_l_integrand,0,a_l,-L_l-d/2,-d/2,'RelTol',int_rtol,'AbsTol',int_atol); 
end

%Upper Cavity End loss (normalised)
if(ANAL_INT)
    P_e_u = real(sum(pi*(a_u/omega/m_0)^2/2.*(abs(A_n.*U_n'.*p_n_u').*besselj(0,h_n_u'.*a_u)).^2,1));
else
    H_rho_u = @(r) sum(-1./(i*omega.*m_0).*p_n_u.'.*A_n.*U_n'.*besselj(1,h_n_u'.*r),1);
    P_e_u = 1/2*2*pi*integral(@(r) r.*abs(H_rho_u(r)).^2,0,a_u,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
end

% Lower Cavity End plate loss (normalised)
if(ANAL_INT)
    P_e_l = real(sum(pi*(a_l/omega/m_0)^2/2.*(abs(D_n.*L_n'.*p_n_l').*besselj(0,h_n_l'.*a_l)).^2,1));
else
    H_rho_l = @(r) sum(1./(i*omega.*m_0).*p_n_l.'.*D_n.*L_n'.*besselj(1,h_n_l'.*r),1);
    P_e_l = 1/2*2*pi*integral(@(r) r.*abs(H_rho_l(r)).^2,0,a_l,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
end

% Upper cavity wall loss (normalised)
if(ANAL_INT)
    A=zeros(N_u,N_u);
    for k = 1:N_u
        for l = 1:N_u
            if(k==l)
                if(isreal(p_n_u(k)))
                    A(k,l)=L_u/2-sin(2*p_n_u(k)*L_u)/(4*p_n_u(k));
                else
                    A(k,l)=-L_u/2+sin(2*p_n_u(k)*L_u)/(4*p_n_u(k));
                end
            else
               A(k,l)=(sin((p_n_u(k)-conj(p_n_u(l)))*L_u)/2/(p_n_u(k)-conj(p_n_u(l))))-(sin((p_n_u(k)+conj(p_n_u(l)))*L_u)/2/(p_n_u(k)+conj(p_n_u(l))));
            end
        end
    end
    P_w_u=real(sum(sum(pi.*a_u/(omega*m_0)^2*(A_n*A_n').*(U_n.'*conj(U_n)).*(h_n_u.'*conj(h_n_u)).*(besselj(0,h_n_u'.*a_u)*besselj(0,h_n_u.*a_u)).*A,1),2));
else
    H_z_u = @(z) sum(-1./(i.*omega.*m_0).*A_n.*U_n.'.*h_n_u.'.*besselj(0,h_n_u.'.*a_u).*sin(p_n_u.'.*(L_u+d/2-z)),1);
    P_w_u = 1./2.*2.*pi.*integral(@(z) a_u.*abs(H_z_u(z)).^2,d/2,L_u+d/2,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
end

% Lower cavity wall loss (normalised)

if(ANAL_INT)
    B=zeros(N_l,N_l);
    for k = 1:N_l
        for l = 1:N_l
            if(k==l)
                if(isreal(p_n_l(k)))
                    B(k,l)=L_l/2-sin(2*p_n_l(k)*L_l)/(4*p_n_l(k));
                else
                    B(k,l)=-L_l/2+sin(2*p_n_l(k)*L_l)/(4*p_n_l(k));
                end
            else
               B(k,l)=(sin((p_n_l(k)-conj(p_n_l(l)))*L_l)/2/(p_n_l(k)-conj(p_n_l(l))))-(sin((p_n_l(k)+conj(p_n_l(l)))*L_l)/2/(p_n_l(k)+conj(p_n_l(l))));
            end
        end
    end
    P_w_l = real(sum(sum(pi.*a_l/(omega*m_0)^2*(D_n*D_n').*(L_n.'*conj(L_n)).*(h_n_l.'*conj(h_n_l)).*(besselj(0,h_n_l'.*a_l)*besselj(0,h_n_l.*a_l)).*B,1),2));
else
    H_z_l = @(z) sum(-1./(i.*omega.*m_0).*D_n.*L_n.'.*h_n_l.'.*besselj(0,h_n_l.'.*a_l).*sin(p_n_l.'.*(L_l+d/2+z)),1);
    P_w_l = 1./2.*2.*pi.*integral(@(z) a_l.*abs(H_z_l(z)).^2,-L_l-d/2,-d/2,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
end

if(ANAL_INT)
    C=zeros(N_u,N_u);
    for k = 1:N_u
        for l = 1:N_u
            if(k==l)
                C(k,l)=a_u^2/2*besselj(0,h_n_u(k)*a_u)^2-a_l^2/2*besselj(0,h_n_u(k)*a_l)^2-a_l^2/2*besselj(1,h_n_u(k)*a_l)^2+a_l./h_n_u(k)*besselj(0,h_n_u(k)*a_l)*besselj(1,h_n_u(k)*a_l);
            else
                C(k,l)=-a_l/(h_n_u(k)^2-h_n_u(l)^2)*(h_n_u(l)*besselj(1,h_n_u(k)*a_l)*besselj(0,h_n_u(l)*a_l)-h_n_u(k)*besselj(0,h_n_u(k)*a_l)*besselj(1,h_n_u(l)*a_l));
            end
        end
    end
    P_f_u=real(sum(sum(pi/(omega*m_0)^2.*((p_n_u.'.*(-A_n.*(U_n.').*(cos(p_n_u*L_u).')))*((p_n_u.'.*(-A_n.*(U_n.').*(cos(p_n_u*L_u).'))))').*C,1),2));
else
    H_rho_u = @(r) sum(-1./(i*omega.*m_0).*p_n_u.'.*A_n.*U_n'.*besselj(1,h_n_u'.*r).*cos(p_n_u.'.*L_u),1);
    P_f_u = 1/2*2*pi*integral(@(r) r*abs(H_rho_u(r)).^2,a_l,a_u,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
end

switch MODE
    case 2 %% SIGMA MODE
        out = omega*m_0/2*(((omega*2*(W_u+W_l))/(Q*(P_e_l+P_e_u+P_w_l+P_w_u+P_f_u))))^-2;
    case 1 %% Q MODE
        out = (omega*2*(W_u+W_l))/(R_s*(P_e_l+P_e_u+P_w_l+P_w_u+P_f_u));
    otherwise
        error('ERROR!!! - Incorrect variable argument');
end
end

 
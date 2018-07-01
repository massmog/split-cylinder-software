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

function [ out ] = Jloss( resonator, coefficients, null_Z, varargin)
%a_l = resonator.a_l;
a_u = resonator.a_u;
%L_l = resonator.L_l;
L_u = resonator.L_u;
f_r = coefficients.f;
e_r = coefficients.e_r;
d = coefficients.d;
sigma = resonator.sigma;

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
            tand = varargin{k+1};
        case 'tand'
            MODE = 2; %% TAND MODE
            Q = varargin{k+1};
        case 'sigma_f'
            MODE = 3; %% CALIBRATE FLANGE LOSS MODE
            Q = varargin{k+1};
        %case 'a_l'
        %    a_l = varargin{k+1};   
        case 'a_u'
            a_u = varargin{k+1};
        %case 'L_l'
        %    L_l = varargin{k+1};
        case 'L_u'
            L_u = varargin{k+1};
        case 'f'
            f_r = varargin{k+1};
        case 'e_r'
            e_r = varargin{k+1};
        case 'd'
            d = varargin{k+1};
        case 'sigma'
            sigma = varargin{k+1};
    end
end

% Coefficients
U_n = coefficients.U_n;
h_n_u = coefficients.h_n_u;
p_n_u = coefficients.p_n_u;

V_m = coefficients.V_m;
%W_m = coefficients.W_m;
h_m_s = coefficients.h_m_s;
p_m_s = coefficients.p_m_s;

%L_n = coefficients.L_n;
%h_n_l =  coefficients.h_n_l;
%p_n_l = coefficients.p_n_l;

N_s = coefficients.N_s;
N_u = coefficients.N_u;
N_l = coefficients.N_u;

% Function tanD - Calculate tangent delta from mode matching problem
b = resonator.b;

A_n = null_Z(1:N_u);
B_m = null_Z(N_u+1:N_u+N_s);
%C_m = null_Z(N_u+N_s+1:N_u+N_s+N_s);
%D_n = null_Z(N_u+N_s+N_s+1:N_u+N_s+N_s+N_u);

omega = 2.*pi.*f_r;
e_0 = 8.85418782e-12;
e_lab = 1.00055;
m_0 = 4*pi*1e-7;
R_s = sqrt(omega*m_0/(2*sigma));
if(MODE==1||MODE==2)
    sigma_f = resonator.sigma_f;
    R_s_f = sqrt(omega*m_0/(2*sigma_f));
end

% Energy in the Sample
E_s = @(r,z) sum(B_m.*V_m.'.*besselj(1,h_m_s.'.*r).*cos(p_m_s.'.*(z)),1);
E_s_integrand = @(r1,z1) arrayfun(@(r,z) r.*abs(E_s(r,z)).^2,r1,z1);

if(ANAL_INT)
    % Integral from z=-d/2 to d/2
    W_s = e_0*e_r*pi*b^2/4.*sum((abs(B_m.*V_m').*besselj(0,h_m_s'.*b)).^2.*(d.*ones(N_s,1)/2+sin(p_m_s.'.*d)./(2*p_m_s.')),1);
else
    W_s = e_0.*e_r./4.*2.*pi.*integral2(E_s_integrand,0,b,-d/2,d/2,'RelTol',int_rtol,'AbsTol',int_atol);
end

% Energy in Resonator
if(ANAL_INT)
    W_u = sum(e_0*e_lab*pi*a_u^2/4.*(abs(A_n.*U_n').*besselj(0,h_n_u'.*a_u)).^2.*((L_u/2)*ones(N_u,1)-sin(2.*p_n_u.'.*L_u)./(4.*p_n_u.')).*(-(imag(p_n_u')~=0)+(imag(p_n_u')==0)),1);
else
    E_u = @(r,z) sum(A_n.*U_n.'.*besselj(1,h_n_u.'.*r).*sin(p_n_u.'.*(L_u+d/2-z)),1);
    E_u_integrand = @(r1,z1) arrayfun(@(r,z) r.*abs(E_u(r,z)).^2,r1,z1);
    W_u = e_0*e_lab./4.*2.*pi.*integral2(E_u_integrand,0,a_u,d/2,L_u+d/2,'RelTol',int_rtol,'AbsTol',int_atol); 
end

%Upper Cavity End plate loss (normalised)
if(ANAL_INT)
    P_e_u = real(sum(pi*(a_u/omega/m_0)^2/2.*(abs(A_n.*U_n'.*p_n_u').*besselj(0,h_n_u'.*a_u)).^2,1));
else
    H_rho_u = @(r) sum(-1./(i*omega.*m_0).*p_n_u.'.*A_n.*U_n'.*besselj(1,h_n_u'.*r),1);
    P_e_u = 1/2*2*pi*integral(@(r) r.*abs(H_rho_u(r)).^2,0,a_u,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
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

% Upper sample region flange loss (normalised)
if(ANAL_INT)
    C=zeros(N_s,N_s);
    for k = 1:N_s
        for l = 1:N_s
            if(k==l)
                C(k,l)=b^2/2*besselj(0,h_m_s(k)*b)^2-a_u^2/2*besselj(0,h_m_s(k)*a_u)^2-a_u^2/2*besselj(1,h_m_s(k)*a_u)^2+a_u./h_m_s(k)*besselj(0,h_m_s(k)*a_u)*besselj(1,h_m_s(k)*a_u);
            else
                C(k,l)=-a_u/(h_m_s(k)^2-h_m_s(l)^2)*(h_m_s(l)*besselj(1,h_m_s(k)*a_u)*besselj(0,h_m_s(l)*a_u)-h_m_s(k)*besselj(0,h_m_s(k)*a_u)*besselj(1,h_m_s(l)*a_u));
            end
        end
    end
    P_f_u=real(sum(sum(pi/(omega*m_0)^2.*((p_m_s.'.*(-B_m.*(V_m.').*(sin(p_m_s*d/2).')))*((p_m_s.'.*(-B_m.*(V_m.').*(sin(p_m_s*d/2).'))))').*C,1),2));
else
    P_f_u = 1/2*2*pi*integral(@(r) r*abs(H_rho_s(r,d/2)).^2,a_u,b,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);
end

% Sample loss
if(ANAL_INT)
    P_s_dash = omega*e_0*e_r*pi*b^2/2.*sum((abs(B_m.*V_m').*besselj(0,h_m_s'.*b)).^2.*(d.*ones(N_s,1)/2+sin(p_m_s.'.*d)./(2*p_m_s.')),1);
else
    P_s_dash = omega.*e_0.*e_r/2*2*pi*integral2(E_s_integrand,0,b,-d/2,d/2,'RelTol',int_rtol,'AbsTol',int_atol);
end
 
% Estimated radiation loss / Poynting vector and power flux
%H_z_s = @(r,z) sum(-h_m_s.'/(1i.*omega.*m_0).*besselj(0,h_m_s.'.*r).*(B_m.*V_m.'.*cos(p_m_s.'*z)+C_m.*W_m.'.*sin(p_m_s.'*z)),1);
%P_rad_integrand = @(r,z) 1/2*E_s(r,z).*conj(H_z_s(r,z));
%P_rad = integral(@(z) 2*pi*b/2*P_rad_integrand(b/2,z),-d/2,d/2,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol);

% Outer-diameter loss
%P_od = 1/2*2*pi*integral(@(z) b*abs(H_z_s(b,z)).^2,-d/2,d/2,'ArrayValued',true,'RelTol',int_rtol,'AbsTol',int_atol)

switch MODE
    case 1 % Q MODE
        out = (omega*2*(W_s+2*W_u))/(R_s*(2*P_e_u+2*P_w_u)+R_s_f*(2*P_f_u)+tand*P_s_dash);
    case 2 % TAND MODE
        out = ((omega*2*(W_s+2*W_u))/Q-R_s*(2*P_e_u+2*P_w_u)-R_s_f*(2*P_f_u))/P_s_dash;
    case 3 % CALIBRATE FLANGE LOSS
        out = omega*m_0/2*((2*P_f_u)/(((omega*2*(W_s+2*W_u))/Q)-R_s*(2*P_e_u+2*P_w_u)))^2;
    otherwise
        error('ERROR!!! - Incorrect variable argument');
end
end

 
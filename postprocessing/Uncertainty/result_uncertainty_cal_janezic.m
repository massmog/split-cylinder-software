function [c_i_table] = result_uncertainty_cal_janezic()

% Includes
addpath('../../scsoft_m12')
addpath('../../lib');
addpath('../../qfactor');
addpath('../..');

% Measurement data
f_r=1.003996067158024e+10;
Q=1.229949696469778e+04;

%% Constants
resonator = constants('janezic'); 

a_u = resonator.a_u;
L_u = resonator.L_u;

%% Find zero
e_r = 1.00055; % Permittivity of air
c_0 = 299792458; % Speed of light [m/s]
m_0 = 4*pi*1e-7; % Permeability of free space


p01_dash=p_nm_dash(0,1);

f_r_a = @(a_0,L_0) c_0/(2*pi*sqrt(e_r))*sqrt((p01_dash./a_0).^2+(1*pi/(2*L_0)).^2);
f_a = @(f_0,L_0) fzero(@(a_0) f_0-f_r_a(a_0, L_0),resonator.a_u);
a = f_a(f_r,L_u);

Q_sigma = @(f_0,a_0,L_0,sigma_0)  c_0/f_0/sqrt(e_r)*(sqrt(2*pi*f_0*m_0*sigma_0/2)).*(p01_dash^2+(1*pi*a_0/(2*L_0))^2)^(3/2)/(2*pi*(p01_dash^2+a_0/L_0*(1*pi*a_0/(2*L_0))^2));
f_sigma = @(f_0,a_0,L_0,Q_0) fzero(@(sigma_0) Q_sigma(f_0,a_0,L_0,sigma_0)-Q_0,3e7);
sigma = f_sigma(f_r,a,L_u,Q);

% Measurement values
f = f_r;
%d = coefficients.d;

% Calculate coefficients
dpf = cdiff(@(y) f_a(y,L_u),f);
dpL = cdiff(@(y) f_a(f,y),L_u);

dqf = cdiff(@(y) f_sigma(y,a,L_u,Q),f);
dqQ = cdiff(@(y) f_sigma(f,a,L_u,y),Q);
dqL = cdiff(@(y) f_sigma(f,a,y,Q),L_u);
dqa = cdiff(@(y) f_sigma(f,y,L_u,Q),a); 

row_names = {'df','dQ','da','dL'};
da_l = [dpf,0,0,dpL]'; 
dsigma  = [dqf,dqQ,dqa,dqL]'; 

c_i_table = table(da_l,dsigma,'RowNames',row_names);
save('c_i_table_cal_janezic.mat','c_i_table');

    % Zeros of the derivative of J_n(x)
    function [ out ] = p_nm_dash( n, m )
        d_bessel = @(n, x) 1./2.*(besselj(n-1,x)-besselj(n+1,x));

        out = zeros(n,m);
        for l = 0:n
            temp = 0;
            for k = 1:100
                 if(sign(d_bessel(l,k*pi/2))~=sign(d_bessel(l,(k+1)*pi/2)))
                     [zero,~,exitflag,~] = fzero(@(x) d_bessel(l,x),[k k+1].*pi/2);
                     if(exitflag)
                        temp = temp + 1;
                        out(l+1,temp)=zero;
                        if(temp>=m)
                            break;
                        end
                     end
                 end
            end
        end

    end

end
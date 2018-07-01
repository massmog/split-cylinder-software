function [] = CAL_Rep_Test_janezic()
%% INCLUDES
addpath('../../scsoft_m12')
addpath('../../lib');
addpath('../../qfactor');
addpath('../..');
%% Constants
resonator = constants('janezic'); 
e_r = 1.00055; % Permittivity of air
c_0 = 299792458; % Speed of light [m/s]
m_0 = 4*pi*1e-7; % Permeability of free space
%% Folder names
testname = '16-Nov-2016-CAL-Rep-Test';
measurement_path = fullfile(fileparts(fileparts(pwd)),'calibrations',testname);
%% Get measurement files
cwd = pwd;
cd(measurement_path);
files = ls('mode*');
cd(cwd);

n_results = size(files,1);

sigma = zeros(1,n_results);
a_l = zeros(1,n_results);
f_r = zeros(1,n_results);
Q = zeros(1,n_results);
f = [];
S21 = [];

for k=1:n_results
    load(fullfile(measurement_path,files(k,:)));
    % Sort by Q factor
    f_Q_S =  sortrows(find_resonances(f,S21,'circlefit'),2);
 
    f_r(k) = f_Q_S(1,1);
    Q(k) = f_Q_S(1,2);

    p01_dash=p_nm_dash(0,1);
    f_r_a = @(a_0,L_0) c_0/(2*pi*sqrt(e_r))*sqrt((p01_dash./a_0).^2+(1*pi/(2*L_0)).^2);
    f_a = @(f_0,L_0) fzero(@(a_0) f_0-f_r_a(a_0, L_0),resonator.a_u);
    
    a_l(k) = f_a(f_Q_S(1,1),resonator.L_u);
    resonator_temp = resonator;
    resonator_temp.a_l = a_l(k);
    
    Q_sigma = @(f_0,a_0,L_0,sigma_0)  c_0/f_0/sqrt(e_r)*(sqrt(2*pi*f_0*m_0*sigma_0/2)).*(p01_dash^2+(1*pi*a_0/(2*L_0))^2)^(3/2)/(2*pi*(p01_dash^2+a_0/L_0*(1*pi*a_0/(2*L_0))^2));
    f_sigma = @(f_0,a_0,L_0,Q_0) fzero(@(sigma_0) Q_sigma(f_0,a_0,L_0,sigma_0)-Q_0,3e7);
    
    sigma(k) = f_sigma(f_Q_S(1,1),a_l(k),resonator.L_u,f_Q_S(1,2));
end
clear('f','S21');

mean(a_l)
sqrt(var(a_l))

mean(sigma)
sqrt(var(sigma))

mean(f_r)
sqrt(var(f_r))

mean(Q)
sqrt(var(Q))

figure
yyaxis left
plot(1:n_results,a_l/a_l(1));
yyaxis right
plot(1:n_results,sigma/sigma(1));

save('results_janezic.mat','a_l','sigma','f_r','Q')
% Zeros of the derivative of J_n(x)
    function [ out ] = p_nm_dash( n, m )
        d_bessel = @(n, x) 1./2.*(besselj(n-1,x)-besselj(n+1,x));

        out = zeros(n,m);
        for l = 0:n
            temp = 0;
            for n = 1:100
                 if(sign(d_bessel(l,n*pi/2))~=sign(d_bessel(l,(n+1)*pi/2)))
                     [zero,~,exitflag,~] = fzero(@(x) d_bessel(l,x),[n n+1].*pi/2);
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

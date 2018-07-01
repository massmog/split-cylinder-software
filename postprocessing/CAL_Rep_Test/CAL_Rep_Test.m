% Prep the workspace
clear
clc
close all
% Includes
addpath('../../scsoft_m12')
addpath('../../lib');
addpath('../../qfactor');
addpath('../..');
%% Constants
resonator = constants(); 
%% Folder names
testname = '16-Nov-2016-CAL-Rep-Test';
measurement_path = fullfile(fileparts(fileparts(pwd)),'calibrations',testname);
%% Get measurement files
cwd = pwd;
cd(measurement_path);
files = ls('mode*');
cd(cwd);

n_results = size(files,1);

sigma = zeros(n_results,1);
a_l = zeros(n_results,1);
f_r = zeros(n_results,1);
Q = zeros(n_results,1);

for k=1:n_results
    load(fullfile(measurement_path,files(k,:)));
    % Sort by Q factor
    f_Q_S =  sortrows(find_resonances(f,S21,'circlefit'),2);
    cal_zer = find_zeros(resonator, @(x) Cmat(resonator, f_Q_S(1,1),'a_l',x),50,[resonator.a_l-1e-4 resonator.a_l+1e-4],'m12','cal','nodisp');
    
    f_r(k) = f_Q_S(1,1);
    Q(k) = f_Q_S(1,2);
    
    a_l(k) = cal_zer{1,1};
    resonator_temp = resonator;
    resonator_temp.a_l = a_l(k);
    
    coefficients = cal_zer{3,1};
    null_Z = cal_zer{4,1};

    sigma(k) = Closs( resonator, coefficients, null_Z, 'sigma',f_Q_S(1,2)); 
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

save('results.mat','a_l','sigma','f_r','Q')

writetable(table((1:n_results)',a_l,sigma,mean(a_l)*ones(n_results,1),mean(sigma)*ones(n_results,1),'VariableNames',{'n' 'al' 'sigma','mual','musigma'}),'reptest_cal.csv')

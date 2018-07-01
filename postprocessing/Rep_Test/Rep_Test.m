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
testname = '16-Nov-2016-PTFE-Rep-Test';
measurement_path = fullfile(fileparts(fileparts(pwd)),'measurements',testname);
%% Get measurement files
cwd = pwd;
cd(measurement_path);
files = ls('mode*');
cd(cwd);
%% Measurement values
d = 1.509E-03;

n_results = size(files,1);

tand = zeros(n_results,1);
e_r = zeros(n_results,1);
f_r = zeros(n_results,1);
Q = zeros(n_results,1);

for k=1:n_results
    load(fullfile(measurement_path,files(k,:)));
    % Sort by Q factor
    f_Q_S =  find_resonances(f,S21,'circlefit');
    [~,r_min] = min(abs(f_Q_S(:,1)-(f(end)+f(1))/2));
    f_r(k) = f_Q_S(r_min,1);
    Q(k) = f_Q_S(r_min,2);
    
    zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_r(k),d),50,[2 2.1]);
    
    
    e_r(k) = zer{1,1};
    resonator_temp = resonator;
    resonator_temp.e_r = e_r(k);
    
    coefficients = zer{3,1};
    null_Z = zer{4,1};

    tand(k) = Zloss( resonator, coefficients, null_Z, 'tand',Q(k)); 
end
clear('f','S21');

mean(e_r)
sqrt(var(e_r))

mean(tand)
sqrt(var(tand))

mean(f_r)
sqrt(var(f_r))

mean(Q)
sqrt(var(Q))

figure
yyaxis left
plot(1:n_results,e_r/e_r(1));
yyaxis right
plot(1:n_results,tand/tand(1));

save('results.mat','e_r','tand','f_r','Q')

writetable(table((1:n_results)',e_r,tand,mean(e_r)*ones(n_results,1),mean(tand)*ones(n_results,1),'VariableNames',{'n' 'er' 'tand','muer','mutand'}),'reptest.csv')


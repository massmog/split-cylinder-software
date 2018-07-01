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
resonator.N = 30;
%% Folder names
testname = '16-Nov-2016-CAL-Broadband';
measurement_path = fullfile(fileparts(fileparts(pwd)),'calibrations',testname);

%% Get measurement files
cwd = pwd;
cd(measurement_path);
files = ls('mode-1*');
cd(cwd);

n_results = size(files,1);

f_max = 20e9;
zer = find_zeros(resonator, @(x) Cmat(resonator, x),1000,[10e9 f_max],'m12','cal','nodisp');
n_zer = size(zer,2);

cal_table = zeros(6,n_zer);

for k=1:n_results
    fprintf('Mode %i:\n',k);
    load(fullfile(measurement_path,files(k,:)));
    % Sort by Q factor
    f_Q_S =  find_resonances(f,S21,'circlefit','finder');
    [~,r_min] = min(abs(f_Q_S(:,1)-f(round(length(f)/2))));
    
    f_r = f_Q_S(r_min,1);
    Q = f_Q_S(r_min,2);
    
    [~,n_min]=min(abs([zer{1,:}]-f_r));
    
    % Plotting
    [~,r_min] = min(abs(f_Q_S(:,1)-(f(end)+f(1))/2));
    [~,x_min] = min(abs(f-f_Q_S(r_min,1))); 
    figure;
    plot(f,20*log10(abs(S21)),'b',f(x_min),20*log10(abs(S21(x_min))),'*');
    text(f(x_min),20*log10(abs(S21(x_min))),sprintf('f_r=%.3e,Q=%.3f',f_Q_S(r_min,1),f_Q_S(r_min,2)),'HorizontalAlignment','right');
   
    
    cal_zer = find_zeros(resonator, @(x) Cmat(resonator, f_r,'a_l',x),50,[resonator.a_l-2e-4 resonator.a_l+2e-4],'m12','cal');
    % Take the mode identical to the mode in question
    % Most likely the zero-finding algorithm suffers from 
    [mode_d,i_min]=min(sum(abs(reshape([cal_zer{2,:}],[2 size(cal_zer,2)])-[zer{2,n_min}]'*ones(1,size(cal_zer,2))),1));
    if(mode_d~=0)
        continue;
    end    
    a_l = cal_zer{1,i_min};
        
    resonator_temp = resonator;
    resonator_temp.a_l = a_l;
    coefficients = cal_zer{3,i_min};
    null_Z = cal_zer{4,i_min};
    mn_mode = cal_zer{2,i_min}';
    

    sigma = Closs( resonator_temp, coefficients, null_Z, 'sigma',Q); 
    cal_table(:,k)=[mn_mode;a_l;sigma;f_r;Q];
end

figure
% a_l Plot
yyaxis left
errorbar(cal_table(5,:),cal_table(3,:),2*1.74e-6*ones(1,size(cal_table,2)))
yyaxis right
errorbar(cal_table(5,:),cal_table(4,:),2*6.3e5*ones(1,size(cal_table,2)))

save('results.mat','cal_table');
mode_names = arrayfun(@(x,y) sprintf('TE%i%i',x,y),cal_table(1,:)',cal_table(2,:)','UniformOutput',false);
writetable(table(mode_names,cal_table(5,:)',cal_table(3,:)',cal_table(4,:)','VariableNames',{'modes','fr','al','sigma'}),'cal_table.csv');
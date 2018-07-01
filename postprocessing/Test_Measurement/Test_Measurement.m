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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Folder names
testname = '17-Nov-2016-HDPE-white';
run_no = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
measurement_path = fullfile(fileparts(fileparts(pwd)),'measurements',testname);

%% Get measurement files
cwd = pwd;
cd(measurement_path);
files = ls(strcat('mode-',sprintf('%i',run_no),'-*'));
cd(cwd);

n_results = size(files,1);

% CALIBRATION TABLE % CORRECTED CALIBRATION!!!!
%cal_table=[1 1 1 1 2 2 2;1 2 3 5 1 2 3;0.0190556261164045 0.0190656744452858 0.0190800124944945 0.0191594068905003 0.0190607688045489 0.0190550560992787 0.0190579595269852;10372062.0682003 9445612.82545431 10144048.0238632 15619668.3968784 9576831.24289870 8165956.58323989 8990930.48346034;10039781631.2110 11298116324.7021 13130450711.1067 17761974604.6582 17796756536.7132 18539867525.9034 19710512418.7658;12486.5485130386 13261.5305753360 15560.1020826924 24026.9976974337 15744.9776807063 15111.5820876677 16773.8201805497];
% Modified TE015 mode
cal_table=[1 1 1 1 2 2 2;1 2 3 5 1 2 3;0.0190556261164045 0.0190656744452858 0.0190800124944945 0.0190800124944945 0.0190607688045489 0.0190550560992787 0.0190579595269852;10372062.0682003 9445612.82545431 10144048.0238632 10144048.0238632 9576831.24289870 8165956.58323989 8990930.48346034;10039781631.2110 11298116324.7021 13130450711.1067 17761974604.6582 17796756536.7132 18539867525.9034 19710512418.7658;12486.5485130386 13261.5305753360 15560.1020826924 24026.9976974337 15744.9776807063 15111.5820876677 16773.8201805497];
% MAXIMUM MEASUREMENT FREQUENCY
f_max = 20e9;

% Load e_r estimate and d from result file
load(fullfile(measurement_path,sprintf('result-%i.mat',run_no))); 
e_r_last = e_r(1);

% Declare variables
e_r = zeros(1,n_results);
tand = zeros(1,n_results);
f_r = zeros(1,n_results);
Q = zeros(1,n_results);
modes = zeros(2,n_results);

for k=1:n_results
    fprintf('Mode %i:\n',k);
    load(fullfile(measurement_path,files(k,:)));
    % Sort by Q factor
    f_Q_S =  find_resonances(f,S21,'circlefit','finder');
    if(isempty(f_Q_S))
        error_msg(sprintf('No resonance curve measured for mode %i - Mode skipped.',k));
        continue;
    end
    % Plotting
    [~,r_min] = min(abs(f_Q_S(:,1)-(f(end)+f(1))/2));
    [~,x_min] = min(abs(f-f_Q_S(r_min,1))); 
    figure;
    plot(f,20*log10(abs(S21)),'b',f(x_min),20*log10(abs(S21(x_min))),'*');
    text(f(x_min),20*log10(abs(S21(x_min))),sprintf('f_r=%.3e,Q=%.3f',f_Q_S(r_min,1),f_Q_S(r_min,2)),'HorizontalAlignment','right');
    % First estimate using general result
    zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_Q_S(r_min,1),d),500,[e_r_last*0.8 1.2*e_r_last]);
    if(isempty(zer))
        continue;
    end
    [~,e_min]=min(abs([zer{1,:}]-e_r_last));
    [mode_d,i_min]=min(sum(abs(cal_table(1:2,:)-[zer{2,e_min}]'*ones(1,length(cal_table))),1));
    if(mode_d~=0)
        error_msg(sprintf('No calibration available for mode %i - Mode skipped.',k));
        continue;
    end
    resonator_temp = resonator;
    resonator_temp.a_l = cal_table(3,i_min);
    resonator_temp.sigma = cal_table(4,i_min);
    resonator_temp.sigma_f = cal_table(4,i_min);
    % Second estimate using calibration
    zer_p = find_zeros(resonator,@(x) Zmat(resonator_temp,x,f_Q_S(r_min,1),d),500,[e_r_last*0.8 1.2*e_r_last]);
    if(isempty(zer_p))
        continue;
    end
    % Check if second calculation yields the same mode
    [mode_p,p_min]=min(sum(abs(reshape([zer_p{2,:}],[2 size(zer_p,2)])-[zer{2,e_min}]'*ones(1,size(zer_p,2))),1));
    if(mode_p~=0)
        error_msg('Mode calibration error!');
        continue;
    end
    % If correct, save e_r value
    e_r(k) = zer_p{1,p_min};
    e_r_last = zer_p{1,p_min};
    
    coeff = zer{3,p_min};
    null_Z = zer{4,p_min};
    tand(k)=Zloss( resonator_temp, coeff, null_Z,'tand',f_Q_S(r_min,2));
    
    % Save resonance curve parameters
    f_r(k)=f_Q_S(r_min,1);
    Q(k)=f_Q_S(r_min,2);
     % Save mode
    modes(:,k)=[zer{2,e_min}]';
    
end

% Save measurement results in MAT file
save(strcat(testname,'.mat'),'e_r','tand','f_r','Q','d','modes');

table_filename = 'e_r_data.csv';
% Read table, if table not available, create table
try e_r_tab = readtable(table_filename,'Delimiter',';');
catch ME
    if (strcmp(ME.identifier,'MATLAB:readtable:OpenFailed'))
        e_r_tab = table;
   end
end



table_name = cell(size(modes,2),1);
table_name{1,1} = testname;
table_empty = cell(size(modes,2),1);
new_e_r_tab = table(table_name,(10*modes(1,:)+modes(2,:))',f_r',e_r',tand',table_empty,'VariableNames',{'Name' 'Mode' 'f_r' 'e_r','tand','comments'});
    
e_r_tab = [e_r_tab;new_e_r_tab];

writetable(e_r_tab,table_filename,'Delimiter',';');

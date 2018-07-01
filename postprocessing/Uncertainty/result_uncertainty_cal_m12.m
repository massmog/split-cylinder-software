function [c_i_table] = result_uncertainty_cal_m12()

% Includes
addpath('../../scsoft_m12')
addpath('../../lib');
addpath('../../qfactor');
addpath('../..');

% Measurement data
f_r=1.003996067158024e+10;
Q=1.229949696469778e+04;

%% Constants
resonator = constants(); 

a_l = resonator.a_l;
a_u = resonator.a_u;
L_u = resonator.L_u;
L_l = resonator.L_l;

%% Find zero
zer_p = find_zeros(resonator, @(x) Cmat(resonator, f_r,'a_l',x),100,[resonator.a_l-1e-4 resonator.a_l+1e-4]);

coefficients = zer_p{3};
null_Z = zer_p{4};


% Measurement values
f = coefficients.f;
d = coefficients.d;

% Calculate coefficients
% DANGER: CORRECT MODES?!
def = cdiff(@(y) find_zeros(resonator, @(x) Cmat(resonator, y,'a_l',x),100,[resonator.a_l-1e-4 resonator.a_l+1e-4],'single','nodisp'),f);
deL_u = cdiff(@(y) find_zeros(resonator, @(x) Cmat(resonator, f,'a_l',x,'L_u',y),100,[resonator.a_l-1e-4 resonator.a_l+1e-4],'single','nodisp'),resonator.L_u);
deL_l = cdiff(@(y) find_zeros(resonator, @(x) Cmat(resonator, f,'a_l',x,'L_l',y),100,[resonator.a_l-1e-4 resonator.a_l+1e-4],'single','nodisp'),resonator.L_l);
dea_u = cdiff(@(y) find_zeros(resonator, @(x) Cmat(resonator, f,'a_l',x,'a_u',y),100,[resonator.a_l-1e-4 resonator.a_l+1e-4],'single','nodisp'),resonator.a_u);

dff = cdiff(@(y) Closs( resonator, coefficients, null_Z, 'sigma',Q,'f',y),f);
dfQ= cdiff(@(y) Closs( resonator, coefficients, null_Z, 'sigma',y),Q);
dfL_u= cdiff(@(y) Closs( resonator, coefficients, null_Z, 'sigma',Q,'L_u',y),resonator.L_u);
dfL_l= cdiff(@(y) Closs( resonator, coefficients, null_Z, 'sigma',Q,'L_l',y),resonator.L_l);
dfa_l = cdiff(@(y) Closs( resonator, coefficients, null_Z, 'sigma',Q,'a_l',y),resonator.a_l); 
dfa_u = cdiff(@(y) Closs( resonator, coefficients, null_Z, 'sigma',Q,'a_u',y),resonator.a_u); 

row_names = {'df','dQ','da_u','da_l','dL_u','dL_l'};
da_l = [def,0,dea_u,0,deL_u,deL_l]'; 
dsigma  = [dff,dfQ,dfa_u,dfa_l,dfL_u,dfL_l]'; 

c_i_table = table(da_l,dsigma,'RowNames',row_names);
save('c_i_table_cal_m12.mat','c_i_table');
end
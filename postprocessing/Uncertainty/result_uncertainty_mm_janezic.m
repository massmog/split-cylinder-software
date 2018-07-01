function [c_i_table] = result_uncertainty_mm_janezic()
% Includes
addpath('../../scsoft_m12')
addpath('../../lib');
addpath('../../qfactor');
addpath('../..');

% Measurement data
f_r=9.661909142153076e+09;
Q=8.881701488892697e+03;
d=1.50915e-3;
e_r_last=2.056259289796823;


%% Constants
resonator = constants(); 

a_l = resonator.a_l;
a_u = resonator.a_u;
L_u = resonator.L_u;
L_l = resonator.L_l;

%% Find zero
zer_p = find_zeros(resonator,@(x) Jmat(resonator,x,f_r,d),500,[e_r_last*0.8 1.2*e_r_last]);

coefficients = zer_p{3};
null_Z = zer_p{4};

% Measurement values
e_r = coefficients.e_r;
f = coefficients.f;
d = coefficients.d;

% Calculate coefficients
% Check modes!!!!!!!!
dgf = cdiff(@(y) find_zeros(resonator, @(x) Jmat(resonator, x, y, d),10,[e_r-1e-2 e_r+1e-2],'single','nodisp'),f);
dga_u = cdiff(@(y) find_zeros(resonator, @(x) Jmat(resonator, x, f, d,'a_u',y),10,[e_r-1e-2 e_r+1e-2],'single','nodisp'),a_u);
%dga_l = cdiff(@(y) find_zeros(resonator, @(x) Zmat(resonator, x, f, d,'a_l',y),10,[e_r-1e-2 e_r+1e-2],'single','nodisp'),a_l);
dgL_u = cdiff(@(y) find_zeros(resonator, @(x) Jmat(resonator, x, f, d,'L_u',y),10,[e_r-1e-2 e_r+1e-2],'single','nodisp'),L_u);
%dgL_l = cdiff(@(y) find_zeros(resonator, @(x) Zmat(resonator, x, f, d,'L_l',y),10,[e_r-1e-2 e_r+1e-2],'single','nodisp'),L_l);
dgd = cdiff(@(y) find_zeros(resonator, @(x) Jmat(resonator, x, f, y),10,[e_r-1e-2 e_r+1e-2],'single','nodisp'),d);


dhf = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q, 'f', y),f);
dhQ = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', y),Q);
dha_u = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q,'a_u',y),a_u);
%dha_l = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q,'a_l',y),a_l);
dhL_u = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q,'L_u',y),L_u);
%dhL_l = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q,'L_l',y),L_l);
dhsigma = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q,'sigma',y),resonator.sigma);
dhd = cdiff(@(y) Jloss( resonator, coefficients, null_Z, 'tand', Q,'d',y),d);

row_names = {'df','dhQ','da_u','dL_u','dd','dsigma'};
de_r  = [dgf,0,dga_u,dgL_u,dgd,0]'; 
dtand  = [dhf,dhQ,dha_u,dhL_u,dhd,dhsigma]'; 

c_i_table = table(de_r,dtand,'RowNames',row_names);
save('c_i_table_mm_janezic.mat','c_i_table');
end
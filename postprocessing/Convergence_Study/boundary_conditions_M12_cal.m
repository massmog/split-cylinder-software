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
resonator.N = 250;

%% Measurement values
f_r = 1.0040e+10;

zer = find_zeros(resonator,@(x) Cmat(resonator,f_r,'a_l',x),100,[resonator.a_l-2e-3 resonator.a_l+2e-3]);
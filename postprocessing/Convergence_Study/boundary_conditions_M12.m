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

%% Measurement values
d = 1.509E-03;
f_r = 9.6619e9;
Q = 8754.3;

zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_r,d),50,[2 2.1]);

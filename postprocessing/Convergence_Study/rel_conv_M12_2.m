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
% Ideal ratio according to Illinski condition N_s = 55 => 1.833

%% Measurement values
d = 1.509E-03;
f_r = 9.6619e9;
Q = 8754.3;

% Constants
omega = 2.*pi.*f_r;
c_0 = 299792458;

% Adjusting mode ratio
N_modes = 70:-1:30;
N_s = min(N_modes):max(N_modes);
e_r = zeros(length(N_modes),1);
tand = zeros(length(N_modes),1);
cond_Z = zeros(length(N_modes),1);

for l = 1:length(N_modes)
resonator.N = N_modes(l);
        for n_iv = 5
            zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_r,d,'N_s',N_s(l)),n_iv,[2.04 2.07],'m12','nodisp');
            if(~isempty(zer))
                break;
            end
        end
        %%%%%%%%%%%%%%%
        if(~isempty(zer))
            coefficients = zer{3,1};
            null_Z = zer{4,1};
            e_r(l) = zer{1,1};
            tand(l) = Zloss( resonator, coefficients, null_Z, 'tand',Q);
            cond_Z(l) = cond(Zmat(resonator,zer{1,1},f_r,d,'N_s',N_s(l)),inf);
        else
            e_r(l) = NaN;
            tand(l) = NaN;
            cond_Z(l) = NaN;
        end
end

plot(N_s./N_modes,10*log10(cond_Z))

%figure
%contour(N_s,N_modes,log10(cond_Z),20)
%hold on
%plot(N_s_ideal,N_modes,'--','color', 'green')

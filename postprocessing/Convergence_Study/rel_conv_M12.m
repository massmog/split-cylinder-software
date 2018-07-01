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
N_modes = 10:1:90;
N_s = 25:round((2*max(N_modes)-25)/length(N_modes)):2*max(N_modes);
e_r = zeros(length(N_modes),length(N_s));
tand = zeros(length(N_modes),length(N_s));
cond_Z = zeros(length(N_modes),length(N_s));
N_s_ideal = zeros(1,length(N_modes));
for l = 1:length(N_modes)
resonator.N = N_modes(l);
    for k = 1:length(N_s)
        for n_iv = 5
            zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_r,d,'N_s',N_s(k)),n_iv,[2.04 2.07],'m12','nodisp');
            if(~isempty(zer))
                break;
            end
        end
        % Progress Sign
        disp(sprintf('%.2f %%',((length(N_s)*(l-1)+k)/(length(N_s)*length(N_modes)))*100))
        %%%%%%%%%%%%%%%
        if(~isempty(zer))
            coefficients = zer{3,1};
            null_Z = zer{4,1};
            e_r(l,k) = zer{1,1};
            tand(l,k) = Zloss( resonator, coefficients, null_Z, 'tand',Q);
            cond_Z(l,k) = cond(Zmat(resonator,zer{1,1},f_r,d,'N_s',N_s(k)));
        else
            e_r(l,k) = NaN;
            tand(l,k) = NaN;
            cond_Z(l,k) = NaN;
        end
    end
    % Compute ideal mode ratio for every N
    [d_p,N_s_ideal(l)]=min(abs(imag((sqrt((omega./c_0.*sqrt(2.056)).^2.-(j1(1:250)./resonator.b).^2)))-imag((sqrt((omega./c_0.*sqrt(1.00055)).^2.-(j1(N_modes(l))./resonator.a_u).^2)))));
    N_s_ideal(l) = N_s_ideal(l)+1*(d_p<0);
end

set(gcf, 'defaultAxesTickLabelInterpreter','latex');
set(gcf, 'defaultLegendInterpreter','latex');

% Plotting
figure
ax=axes;
e_r_c = min(min(e_r)):(max(max(e_r))-min(min(e_r)))/20:max(max(e_r));
[C,h] = contour(N_modes,N_s,e_r',e_r_c);
clabel(C,h,'Interpreter','latex');
clabel(C,h,'FontSize',11);
h.TextList = e_r_c(e_r_c>2.053);
h.TextListMode = 'manual';
h.ShowText = 'on';
hold on
plot(N_modes,N_s_ideal,'--','color', 'red')
set(ax, 'FontSize', 12)
xlabel('$N_u$','fontsize',14)
ylabel('$N_s$','fontsize',14);
print('RC_er.eps','-depsc2')

figure
ax=axes;
[C,h]=contour(N_modes,N_s,tand',15,'ShowText','on')
clabel(C,h,'Interpreter','latex');
clabel(C,h,'FontSize',11);
hold on
plot(N_modes,N_s_ideal,'--','color', 'red')

set(ax, 'FontSize', 12)
xlabel('$N_u$','fontsize',14)
ylabel('$N_s$','fontsize',14);
print('RC_tand.eps','-depsc2');


%figure
%contour(N_s,N_modes,log10(cond_Z),20)
%hold on
%plot(N_s_ideal,N_modes,'--','color', 'green')

save('results_RC.mat');
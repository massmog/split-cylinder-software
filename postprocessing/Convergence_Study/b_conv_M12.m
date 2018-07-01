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
N_modes = [20,30,40,50,60,70,80];
% Ideal ratio according to Illinski condition N_s = 55 => 1.833

%% Measurement values
d = 1.509E-03;
f_r = 9.6619e9;
Q = 8754.3;

% Adjusting mode ratio
b_vector = (resonator.a_u+1e-4):2e-4:45e-3;
e_r = zeros(length(N_modes),length(b_vector));
tand = zeros(length(N_modes),length(b_vector));

for l = 1:length(N_modes)
    resonator.N = N_modes(l);
    for k = 1:length(b_vector)
        resonator.b = b_vector(k);
        zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_r,d),50,[2.0 2.1],'nodisp');
        coefficients = zer{3,1};
        null_Z = zer{4,1};
        if(zer{2,1}~=[1 1])
            break;
        end

        e_r(l,k) = zer{1,1};
        tand(l,k) = Zloss( resonator, coefficients, null_Z, 'tand',Q); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
%mean_e_r = mean(e_r(b_vector>35e-3));
%mean_tand = mean(tand(b_vector>35e-3));
%sigma_e_r = 2*sqrt(var(e_r((b_vector>30e-3)&(b_vector<35e-3))));
%sigma_tand = 2*sqrt(var(tand((b_vector>30e-3)&(b_vector<35e-3))));
%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf, 'defaultAxesTickLabelInterpreter','latex');
set(gcf, 'defaultLegendInterpreter','latex');
ax = axes;
plot(b_vector,e_r(1,:))
hold on
for m = 2:size(e_r,1)
plot(b_vector,e_r(m,:))
%hline = refline([0 mean_e_r]);
end
y_bnd = ylim;
ylim manual
ln1=line([resonator.a_u,resonator.a_u],y_bnd)
set(ln1,'Color','red','LineStyle','--')
ln2=line([35e-3,35e-3],y_bnd)
set(ln2,'Color','red','LineStyle','--')
legend('N=20','N=30','N=40','N=50','N=60','N=70','N=80')
set(ax, 'FontSize', 12)
xlabel('b [m]','Interpreter','latex','fontsize',14)
ylabel('$\epsilon_r$ [1]','Interpreter','latex','fontsize',14);
x_tks = 15e-3:5e-3:45e-3;
x_lbls=num2cell(x_tks);
[~,x_tks_m]=min(abs(resonator.a_u-x_tks))
x_tks(x_tks_m)=resonator.a_u;
x_lbls{x_tks_m}='$a_u$';
[~,x_tks_m]=min(abs(35e-3-x_tks))
x_tks(x_tks_m)=35e-3;
x_lbls{x_tks_m}='$b_{flange}$';
set(ax,'Xtick',x_tks);
set(ax,'XtickLabel',x_lbls);
print('b_conv_er.eps','-depsc2')

%%%%%%%%%%%%%%%%%%%%%%%
figure
set(gcf, 'defaultAxesTickLabelInterpreter','latex');
set(gcf, 'defaultLegendInterpreter','latex');
ax = axes;
plot(b_vector,tand(1,:))
hold on
for m = 2:size(tand,1)
plot(b_vector,tand(m,:))
%hline = refline([0 mean_e_r]);
end
y_bnd = ylim;
ylim manual
ln1=line([resonator.a_u,resonator.a_u],y_bnd)
set(ln1,'Color','red','LineStyle','--')
ln2=line([35e-3,35e-3],y_bnd)
set(ln2,'Color','red','LineStyle','--')
legend('N=20','N=30','N=40','N=50','N=60','N=70','N=80')
set(ax, 'FontSize', 12)
xlabel('b [m]','Interpreter','latex','fontsize',14)
ylabel('$\tan\delta$ [1]','Interpreter','latex','fontsize',14);
x_tks = 15e-3:5e-3:45e-3;
x_lbls=num2cell(x_tks);
[~,x_tks_m]=min(abs(resonator.a_u-x_tks))
x_tks(x_tks_m)=resonator.a_u;
x_lbls{x_tks_m}='$a_u$';
[~,x_tks_m]=min(abs(35e-3-x_tks))
x_tks(x_tks_m)=35e-3;
x_lbls{x_tks_m}='$b_{flange}$';
set(ax,'Xtick',x_tks);
set(ax,'XtickLabel',x_lbls);
print('b_conv_tand.eps','-depsc2')

save('results_b.mat');
% Copyright (C) 2017  Max Gattringer

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

function [E_phi_max, H_z_u_loop, n_mode, E_phi_specimen ] = plot_mode(resonator, coefficients, null_Z, mode_name, varargin)
%% MODEL PARSER
JAN = 0;
CAL = 0;
M12 = 0;
NODISP = 0;

nVarargs = length(varargin);
for k = 1:nVarargs
    switch varargin{k}
        case 'janezic'
            JAN = 1;
        case 'm12'
            M12 = 1;
        case 'cal'
            CAL = 1;
        case 'nodisp'
            NODISP = 1;
    end
end

f = coefficients.f;
d = coefficients.d;
if(~CAL)
    e_r = coefficients.e_r;
end

omega = 2.*pi.*f;
e_0 = 8.85418782e-12;
e_lab = 1.00055;
m_0 = 4*pi*1e-7;
h_feed = 12.5e-3;

U_n = coefficients.U_n;
h_n_u = coefficients.h_n_u;
p_n_u = coefficients.p_n_u;

N_u = coefficients.N_u;

    % JANEZIC
    if(JAN)
        N_s = coefficients.N_s;
        V_m = coefficients.V_m;
        h_m_s = coefficients.h_m_s;
        p_m_s = coefficients.p_m_s;
    end
    % M12
    if(M12&&~CAL)
        N_s = coefficients.N_s;
        V_m = coefficients.V_m;
        h_m_s = coefficients.h_m_s;
        p_m_s = coefficients.p_m_s;
        h_n_l =  coefficients.h_n_l;
        p_n_l = coefficients.p_n_l;
        L_n = coefficients.L_n;
        W_m = coefficients.W_m;
    end
    % M12 CAL
    if(M12&&CAL)
        N_s = 0;
        h_n_l =  coefficients.h_n_l;
        p_n_l = coefficients.p_n_l;
        L_n = coefficients.L_n;
    end



% Plot resonator modes for rho=a/2
L_u = resonator.L_u;
a_u = resonator.a_u;
L_l = resonator.L_l;
a_l = resonator.a_l;
b   = resonator.b;

N_points = 2000;
d_s = 10; % Density of points in substrate region (faster changes!)
N_s_points = round(N_points*d_s*d/(L_l+d_s*d+L_u));
N_u_points = round(N_points*L_l/(L_l+d_s*d+L_u));
N_l_points = round(N_points*L_u/(L_l+d_s*d+L_u));

z_u = d/2:L_u/(N_u_points-1):L_u+d/2;
z_s = -d/2:d/(N_s_points-1):d/2;
z_l = -L_l-d/2:L_l/(N_l_points-1):-d/2;
rho = min([a_u,a_l])-1e-3;

A_n = null_Z(1:N_u);
B_m = null_Z(N_u+1:N_u+N_s);
if(M12)
    C_m = null_Z(N_u+N_s+1:N_u+N_s+N_s);
    D_n = null_Z(N_u+N_s+N_s+1:N_u+N_s+N_s+N_u);
end

E_phi_u = sum((A_n.*U_n.'.*besselj(1,h_n_u*rho)')*ones(1,length(z_u)).*sin(p_n_u.'*(L_u+d/2-z_u)),1);
if(~CAL)
    if(M12)
        E_phi_s = sum((B_m.*V_m.'.*besselj(1,h_m_s*rho)')*ones(1,length(z_s)).*cos(p_m_s.'*(z_s)),1)+sum((C_m.*W_m.'.*besselj(1,h_m_s*rho)')*ones(1,length(z_s)).*sin(p_m_s.'*(z_s)),1);
    end
    if(JAN)
        E_phi_s = sum((B_m.*V_m.'.*besselj(1,h_m_s*rho)')*ones(1,length(z_s)).*cos(p_m_s.'*(z_s)),1);
    end
else
    E_phi_s = [];
end

if(M12)
    E_phi_l = sum((D_n.*L_n.'.*besselj(1,h_n_l*rho)')*ones(1,length(z_l)).*sin(p_n_l.'*(L_l+d/2+z_l)),1);
end

%% Normalization constants
% M12
if(M12&&~CAL)
    V=pi*b^2*d+pi*a_l^2*L_l+pi*a_u^2*L_u;
    % Energy in the Sample
    W_s = e_0*e_r*pi*b^2/4.*sum((abs(B_m.*V_m').*besselj(0,h_m_s'.*b)).^2.*(d.*ones(N_s,1)/2+sin(p_m_s.'.*d)./(2*p_m_s.'))+(-(imag(p_m_s')~=0)+(imag(p_m_s')==0)).*(abs(C_m.*W_m').*besselj(0,h_m_s'.*b)).^2.*(d.*ones(N_s,1)/2-sin(p_m_s.'.*d)./(2*p_m_s.')),1);
    % Energy in the Upper Cavity
    W_u = sum(e_0*e_lab*pi*a_u^2/4.*(abs(A_n.*U_n').*besselj(0,h_n_u'.*a_u)).^2.*((L_u/2)*ones(N_u,1)-sin(2.*p_n_u.'.*L_u)./(4.*p_n_u.')).*(-(imag(p_n_u')~=0)+(imag(p_n_u')==0)),1);
    % Energy in the Lower Cavity
    W_l = sum(e_0*e_lab*pi*a_l^2/4.*(abs(D_n.*L_n').*besselj(0,h_n_l'.*a_l)).^2.*((L_l/2)*ones(N_u,1)-sin(2.*p_n_l.'.*L_l)./(4.*p_n_l.')).*(-(imag(p_n_l')~=0)+(imag(p_n_l')==0)),1);
    W_e = W_s+W_u+W_l;    
    k_e = sqrt((4/e_0/e_lab*(W_u+W_l)+4/e_0/e_r*W_s)/V);
end
if(JAN)
    V=pi*b^2*d+2*pi*a_u^2*L_u;
    % Energy in the sample
    W_s = e_0*e_r*pi*b^2/4.*sum((abs(B_m.*V_m').*besselj(0,h_m_s'.*b)).^2.*(d.*ones(N_s,1)/2+sin(p_m_s.'.*d)./(2*p_m_s.')),1);
    % Energy in the Upper Cavity
    W_u = sum(e_0*e_lab*pi*a_u^2/4.*(abs(A_n.*U_n').*besselj(0,h_n_u'.*a_u)).^2.*((L_u/2)*ones(N_u,1)-sin(2.*p_n_u.'.*L_u)./(4.*p_n_u.')).*(-(imag(p_n_u')~=0)+(imag(p_n_u')==0)),1);
    W_e = W_s+2*W_u;
    k_e = sqrt((4/e_0/e_lab*(2*W_u)+4/e_0/e_r*W_s)/V);
end
if(M12&&CAL)
    V=pi*a_l^2*L_l+pi*a_u^2*L_u;
    % Energy in the Upper Cavity
    W_u = sum(e_0*e_lab*pi*a_u^2/4.*(abs(A_n.*U_n').*besselj(0,h_n_u'.*a_u)).^2.*((L_u/2)*ones(N_u,1)-sin(2.*p_n_u.'.*L_u)./(4.*p_n_u.')).*(-(imag(p_n_u')~=0)+(imag(p_n_u')==0)),1);
    % Energy in the Lower Cavity
    W_l = sum(e_0*e_lab*pi*a_l^2/4.*(abs(D_n.*L_n').*besselj(0,h_n_l'.*a_l)).^2.*((L_l/2)*ones(N_u,1)-sin(2.*p_n_l.'.*L_l)./(4.*p_n_l.')).*(-(imag(p_n_l')~=0)+(imag(p_n_l')==0)),1);
    W_e = W_u+W_l; 
    k_e = sqrt((4/e_0/e_lab*(W_u+W_l))/V);
end
k_m=sqrt(4*W_e/V/m_0);
%% DEBUG FUNCTION - PLOT BOUNDARY VALUES
if(0)
    d_rho= 0.00001;
    rho_b = 0:d_rho:(b-1e-5);
    rho_a_u = 0:d_rho:a_u;
    rho_a_l = 0:d_rho:a_l;
    % E-FIELDS
    E_phi_u_BC = sum((A_n.*U_n.'.*sin(p_n_u*(L_u))'*ones(1,size(rho_a_u,2))).*besselj(1,h_n_u'*rho_a_u),1)./k_e;
    if(~CAL)
    E_phi_s_BC_u = (sum((B_m.*V_m.'.*cos(p_m_s*d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1)+sum((C_m.*W_m.'.*sin(p_m_s*d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1))./k_e;
    E_phi_s_BC_l = (sum((B_m.*V_m.'.*cos(p_m_s*-d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1)+sum((C_m.*W_m.'.*sin(p_m_s*-d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1))./k_e;
    end
    E_phi_l_BC = sum((D_n.*L_n.'.*sin(p_n_l*(L_l))'*ones(1,size(rho_a_l,2))).*besselj(1,h_n_l'*rho_a_l),1)./k_e;
    % H-FIELDS
    H_rho_u_BC = sum((-(p_n_u./(1i*omega*m_0)).'.*A_n.*U_n.'.*cos(p_n_u*(L_u))'*ones(1,size(rho_a_u,2))).*besselj(1,h_n_u'*rho_a_u),1)./k_m;
    if(~CAL)
    H_rho_s_BC_u = (sum((-(p_m_s./(1i*omega*m_0)).'.*B_m.*V_m.'.*sin(p_m_s*d/2)'*ones(1,size(rho_a_u,2))).*besselj(1,h_m_s'*rho_a_u),1)+sum(((p_m_s./(1i*omega*m_0)).'.*C_m.*W_m.'.*cos(p_m_s*d/2)'*ones(1,size(rho_a_u,2))).*besselj(1,h_m_s'*rho_a_u),1))./k_m;
    H_rho_s_BC_l = (sum((-(p_m_s./(1i*omega*m_0)).'.*B_m.*V_m.'.*sin(p_m_s*-d/2)'*ones(1,size(rho_a_l,2))).*besselj(1,h_m_s'*rho_a_l),1)+sum(((p_m_s./(1i*omega*m_0)).'.*C_m.*W_m.'.*cos(p_m_s*-d/2)'*ones(1,size(rho_a_l,2))).*besselj(1,h_m_s'*rho_a_l),1))./k_m;
    %H_rho_s_BC_u = (sum((-(p_m_s./(1i*omega*m_0)).'.*B_m.*V_m.'.*sin(p_m_s*d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1)+sum(((p_m_s./(1i*omega*m_0)).'.*C_m.*W_m.'.*cos(p_m_s*d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1))./k_m;
    %H_rho_s_BC_l = (sum((-(p_m_s./(1i*omega*m_0)).'.*B_m.*V_m.'.*sin(p_m_s*-d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1)+sum(((p_m_s./(1i*omega*m_0)).'.*C_m.*W_m.'.*cos(p_m_s*-d/2)'*ones(1,size(rho_b,2))).*besselj(1,h_m_s'*rho_b),1))./k_m;
    end
    H_rho_l_BC = sum((-(p_n_u./(1i*omega*m_0)).'.*D_n.*L_n.'.*cos(p_n_l*(L_l))'*ones(1,size(rho_a_l,2))).*besselj(1,h_n_l'*rho_a_l),1)./k_m;
    % PLOTS
    if(~CAL)
    figure;
    ax = axes;
    set(gcf, 'defaultAxesTickLabelInterpreter','latex');
    set(gcf, 'defaultLegendInterpreter','latex');
    semilogy(rho_b,(abs(E_phi_s_BC_u)));
    hold on;
    semilogy(rho_a_u,(abs(E_phi_u_BC)));
     ln2=line([resonator.a_u,resonator.a_u],ylim);
    set(ln2,'Color','red','LineStyle','--')
    set(ax, 'FontSize', 12)
    x_tks = 0:5e-3:resonator.b;
    [~,x_tks_m]=min(abs(resonator.a_u-x_tks));
    x_tks(x_tks_m)=resonator.a_u;
    x_lbls=num2cell(x_tks);    
    x_lbls{end}='$b$';
    x_lbls{x_tks_m}='$a_u$';
    set(ax,'Xtick',x_tks);
    set(ax,'XtickLabel',x_lbls);
    xlabel('$\rho$ [m]','Interpreter','latex','fontsize',14)
    ylabel('$|E_{\phi}(\rho,z=\frac{d}{2})|$ [1]','Interpreter','latex','fontsize',14);
    hold off;
    legend('Sample region','Upper cavity')
    print(fullfile(pwd,'E_phi_upper.eps'),'-depsc2');
    
    figure;
    ax = axes;
    set(gcf, 'defaultAxesTickLabelInterpreter','latex');
    set(gcf, 'defaultLegendInterpreter','latex');
    semilogy(rho_b,(abs(E_phi_s_BC_l)));
    hold on;
    semilogy(rho_a_l,(abs(E_phi_l_BC)));
    ln1=line([resonator.a_l,resonator.a_l],ylim);
    set(ln1,'Color','red','LineStyle','--')
    set(ax, 'FontSize', 12)
    x_tks = 0:5e-3:resonator.b;
    [~,x_tks_m]=min(abs(resonator.a_u-x_tks));
    x_tks(x_tks_m)=resonator.a_l;
    x_lbls=num2cell(x_tks);  
    x_lbls{end}='$b$';
    x_lbls{x_tks_m}='$a_l$';
    set(ax,'Xtick',x_tks);
    set(ax,'XtickLabel',x_lbls);
    xlabel('$\rho$ [m]','Interpreter','latex','fontsize',14)
    ylabel('$|E_{\phi}(\rho,z=-\frac{d}{2})|$ [1]','Interpreter','latex','fontsize',14);
    legend('Sample region','Lower cavity')
    print(fullfile(pwd,'E_phi_lower.eps'),'-depsc2');
    
    figure;
    ax = axes;
    set(gcf, 'defaultAxesTickLabelInterpreter','latex');
    set(gcf, 'defaultLegendInterpreter','latex');
    semilogy(rho_a_u,(abs(H_rho_s_BC_u)));
    hold on;
    semilogy(rho_a_u,(abs(H_rho_u_BC)));
    set(ax, 'FontSize', 12)
    x_tks = [0:5e-3:15e-3,resonator.a_u];
    x_lbls=num2cell(x_tks);  
    x_lbls{end}='$a_u$';
    set(ax,'Xtick',x_tks);
    set(ax,'XtickLabel',x_lbls)
    xlabel('$\rho$ [m]','Interpreter','latex','fontsize',14)
    ylabel('$|H_{\rho}(\rho,z=\frac{d}{2})|$ [1]','Interpreter','latex','fontsize',14);
    legend('Sample region','Upper cavity')
    print(fullfile(pwd,'H_rho_upper.eps'),'-depsc2');
    
    figure;
    ax = axes;
    set(gcf, 'defaultAxesTickLabelInterpreter','latex');
    set(gcf, 'defaultLegendInterpreter','latex');
    semilogy(rho_a_u,(abs(H_rho_s_BC_l)));
    hold on;
    semilogy(rho_a_l,(abs(H_rho_l_BC)));
    set(ax, 'FontSize', 12)
    x_tks = [0:5e-3:15e-3,resonator.a_l];
    x_lbls=num2cell(x_tks);  
    x_lbls{end}='$a_l$';
    set(ax,'Xtick',x_tks);
    set(ax,'XtickLabel',x_lbls)
    xlabel('$\rho$ [m]','Interpreter','latex','fontsize',14)
    ylabel('$|H_{\rho}(\rho,z=-\frac{d}{2})|$ [1]','Interpreter','latex','fontsize',14);   
    legend('Sample region','Lower cavity')
    print(fullfile(pwd,'H_rho_lower.eps'),'-depsc2');
    end
    
    if(CAL)
    figure;
    ax = axes;
    set(gcf, 'defaultAxesTickLabelInterpreter','latex');
    set(gcf, 'defaultLegendInterpreter','latex');
    semilogy(rho_a_u,(abs(E_phi_u_BC)));
    hold on;
    semilogy(rho_a_l,(abs(E_phi_l_BC)));
    set(ax, 'FontSize', 12)
    x_tks = [0:5e-3:15e-3,resonator.a_l];
    x_lbls=num2cell(x_tks);  
    x_lbls{end}='$a_l$';
    set(ax,'Xtick',x_tks);
    set(ax,'XtickLabel',x_lbls)
    xlabel('$\rho$ [m]','Interpreter','latex','fontsize',14)
    ylabel('$|E_{\phi}(\rho,z=0)|$ [1]','Interpreter','latex','fontsize',14);   
    legend('Upper cavity','Lower cavity')
    print(fullfile(pwd,'E_phi_cal.eps'),'-depsc2');
    
    figure;
    ax = axes;
    set(gcf, 'defaultAxesTickLabelInterpreter','latex');
    set(gcf, 'defaultLegendInterpreter','latex');
    semilogy(rho_a_l,(abs(H_rho_u_BC(1:length(rho_a_l)))));
    hold on;
    semilogy(rho_a_l,(abs(H_rho_l_BC)));
    set(ax, 'FontSize', 12)
    x_tks = [0:5e-3:15e-3,resonator.a_l];
    x_lbls=num2cell(x_tks);  
    x_lbls{end}='$a_l$';
    set(ax,'Xtick',x_tks);
    set(ax,'XtickLabel',x_lbls)
    xlabel('$\rho$ [m]','Interpreter','latex','fontsize',14)
    ylabel('$|H_{\rho}(\rho,z=0)|$ [1]','Interpreter','latex','fontsize',14);   
    legend('Upper Cavity','Lower cavity')
    print(fullfile(pwd,'H_rho_cal.eps'),'-depsc2');
    end
end

%% Magnetic field strength at coupling loop
h_loop = d/2+h_feed;
H_z_u_loop = sum(-1./(i.*omega.*m_0).*A_n.*U_n.'.*h_n_u.'.*besselj(0,h_n_u.'.*a_u).*sin(p_n_u.'.*(L_u+d/2-h_loop)),1)/k_m;


if(~CAL)
    if(M12)
        z = [z_l z_s z_u];
        E_phi = [E_phi_l E_phi_s E_phi_u]/k_e;
    end
    if(JAN)
        z = [-fliplr(z_u) z_s z_u];
        E_phi = [fliplr(E_phi_u) E_phi_s E_phi_u]/k_e;
    end
else
    z = [z_l z_u];
 	E_phi = [E_phi_l E_phi_u]/k_e;  
end

E_phi_max = max(E_phi);
if(CAL==0)
    E_phi_specimen = E_phi(round(N_points/2));
end

% Estimate mode number from phase zero crossings
n_mode = sum(diff(wrapTo360(round(180/pi*(angle(E_phi(50:end-50))))))~=0)+1;

if(~NODISP)
    figure
    ax = axes;
    %yyaxis left
    set(ax, 'FontSize', 14)
    p1=plot(z, abs(E_phi));
    set(ax,'Xtick',[-L_l-d/2 -h_loop -d/2-3e-3 0 d/2+2e-3 h_loop L_u+d/2])
    set(ax,'XtickLabel',{'$-L_l-d/2$','$-h_l$','$-d/2$','0','$d/2$','$h_l$','$L_u+d/2$'})
    set(ax,'TickLabelInterpreter', 'latex');
    xlabel('z [m]','Interpreter','latex','fontsize',14)
    ylabel('$|E_\phi |$ [1]','Interpreter','latex','fontsize',14);
    yyaxis right
    p2=plot(z, 180/pi*(angle(E_phi)));
    xlabel('$z$','Interpreter','latex','fontsize',14)
    ylabel('$\arg(E_{\phi})$ [o]','Interpreter','latex','fontsize',14);
    ylim([-180 180])   
    set(p1,'linewidth',1.5);
    set(p2,'linewidth',1.5);
    ax.YAxis(2).TickValues=[-180,-90,0,90,180];
    
    leg=legend('$|E_\phi |$','arg$(E_{\phi})$'); 
    set(leg,'Interpreter','latex','fontsize',12);
    ax.XAxis.TickLength = [0,0];
    ax.XAxis.FontSize = 11;
    ax.YAxis(1).FontSize = 11;
    ax.YAxis(2).FontSize = 11;
    % Wall markers
    line([d/2+L_u d/2+L_u],get(ax,'YLim'),'Color',[1 0 0]);
    line([-d/2-L_l -d/2-L_l],get(ax,'YLim'),'Color',[1 0 0]);
    % Substrate region marker
    line([-d/2 -d/2],get(ax,'YLim'),'Color',[1 0 0]);
    line([d/2 d/2],get(ax,'YLim'),'Color',[1 0 0]);

    % Feed marker
    line([d/2+h_feed d/2+h_feed],get(ax,'YLim'),'Color',[1 .5 0],'LineStyle','--');
    line([-d/2-h_feed -d/2-h_feed],get(ax,'YLim'),'Color',[1 .5 0],'LineStyle','--');

    title(sprintf('%s,%i at %.3f GHz, $K_{L}=%.3f$ dB',mode_name,n_mode,coefficients.f/1e9,20*log10(abs(H_z_u_loop))),'Interpreter','latex','fontsize',13);
end

end


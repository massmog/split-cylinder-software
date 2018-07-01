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

function [ f_r, Q] = CircleFit( f_in, S_in, varargin)
    plot_fit = 0;

    nVarargs = length(varargin);
    for k = 1:nVarargs
        switch varargin{k}
           case 'plot_fit'
                plot_fit  = 1;
        end
    end

    % Check input
    if(~isvector(f_in) && ~isvector(S_in))
        error('Input variables must be vectors!');
    end
    % Transform row vectors to column vectors
    if(~iscolumn(f_in))
        f_in    = f_in.';
    end
    
    if(~iscolumn(S_in))
        S_in    = S_in.';        
    end
      
    % Reshape resonance curve for correct fit
    % Frequency spacing factor of Delta = 2.6 (2.6 times resonance bandwidth)
    Delta       = 1.7;
    [f, S]      = ReshapeResonanceCurve(f_in, S_in, Delta);
    % Initial guess
    [f_r_0, Q_0]= ThreeDBmethod( f_in, S_in );
    % Specifying our tranmission coefficient estimator
    S21_est     = @(f,f_r,Q,gamma_1,gamma_2) gamma_1./(1+j.*Q.*(f./f_r-f_r./f))+gamma_2;
    % Define an objective function
    objfcn      = @(v) S21_est(f,v(1),v(2),v(3)+i*v(4),v(5)+i*v(6)) - S;
    % Fitting options    
    opts        = optimoptions(    @lsqnonlin,...
                            'Algorithm','levenberg-marquardt',...
                            'Display','off',...
                            'StepTolerance',1e-50,...
                            'FunctionTolerance',1e-50,...
                            'MaxFunctionEvaluations',1500,...
                            'MaxIterations',500);
    % Initial guess for the fitting
    x_0         = [f_r_0,Q_0,real(max(S_in)),imag(max(S_in)),real(min(S_in)),imag(min(S_in))];
    % Non-linear least square fitting of real and imaginary part of the
    % function
    [x_est,resnorm,residuals,exitflag,output] = lsqnonlin(@(v) [real(objfcn(v)) imag(objfcn(v))],x_0,[],[],opts);
    
    %<DEBUG>
    %figure
    %plot(f,abs(objfcn(x_est)))
    %hold
    %plot(f,abs(objfcn(x_0)))
    
    if(plot_fit)
        gamma_1 = x_est(3)+i*x_est(4);
        gamma_2 = x_est(5)+i*x_est(6);
        figure
        subplot(2,1,1);
        plot(f/1e9,20*log10(abs(S21_est(f,x_est(1),x_est(2),gamma_1,gamma_2))),'r',f/1e9,20*log10(abs(S)),'b');
        subplot(2,1,2);
        plot(f/1e9,real(objfcn(x_est)),'g',f/1e9,imag(objfcn(x_est)),'b');
        %% SNR according to Petersan, Anlage - Measurement of resonant frequency and quality factor of microwave resonators: Comparison of methods
        d_i = abs(S-(gamma_2+0.5*gamma_1));
        r_circle = (abs(gamma_1)/2);
        SNR = r_circle/sqrt(sum((d_i-r_circle).^2,1)/length(d_i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
        title(sprintf('R_{norm} = %.2e, SNR_{est} = %.2f',resnorm, SNR)); 
    end
    %</DEBUG>
    
    f_r         = x_est(1);
    Q           = x_est(2);
end




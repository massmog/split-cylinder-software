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

function [ f_r, Q ] = Coakley( f_in, S_in  )
    % Check input
    if(~isvector(f_in)&&~isvector(S_in))
        error('Input variables must be vectors!');
    end
    % Transform row vectors to column vectors
    if(~iscolumn(f_in))
        f_in = f_in.';
    end
    
    if(~iscolumn(S_in))
        S_in = S_in.';        
    end
    
    % Q-Factor estimation to Coakley et al. (IEEE Trans.MTT, March 2003)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: Compute Q factor using the Estin method:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ f_r_estin, Q_estin, T_estin ] = Estin(f_in, S_in);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2: Use Estin-Estimate as start point for non-linear LS estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reshape resonance curve for correct LS linear fit
    % Frequency spacing factor of Delta = 2.6 (Not optimized!)
    Delta=2.9;
    [f,S] = ReshapeResonanceCurve(f_in, S_in, Delta);
    T=(abs(S)./max(abs(S))).^2;
    % NLLS Estimator
    setupOptions();
    theta_options.StartPoint = [T_estin, Q_estin, f_r_estin, min(T)];

    [NLLS_fit, gof, stats] = fit(f,T,theta_estimator, theta_options);
    if(sqrt(gof.rmse)/norm(T)>0.01)
        warning('First fit failed!');
    end
    NLLS_coeff = coeffvalues(NLLS_fit);
    % figure
    % plot(NLLS_fit, f, T)
    Q_NLLS = NLLS_coeff(2);
    f_r_NLLS = NLLS_coeff(3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 3: Estimate variance function and weights based on binned squared
    % residuals
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = length(f);
    K = round(N/50);
    % Degree of freedom factor (dubious?!)
    DFF = N/(N-4);
    % Refactor length
    r_length = N-rem(N,K);
    % Binning
    residuals = (stats.residuals).^2;
    binned_f = mean(reshape(f(1:r_length),K,r_length/K),1).';
    binned_residuals = DFF.*mean(reshape(residuals(1:r_length),K,r_length/K),1).';
    % Variance function fit
    gamma_estimator = fittype(sprintf('gamma_1/(1+%1$d^2*(x/%2$d-%2$d/x)^2)+gamma_2',Q_NLLS,f_r_NLLS));
    gamma_options = fitoptions(gamma_estimator);
    gamma_options.Lower = [min(binned_residuals)*0.01 min(binned_residuals)*0.01];
    gamma_options.Upper = [max(binned_residuals)*1.5 max(binned_residuals)];
    gamma_options.TolFun = min(binned_residuals)*1e-10;
    gamma_options.TolX = 1e-12; 
    gamma_options.StartPoint = [max(binned_residuals)*0.5 min(binned_residuals)];
    [gamma_fit, ggof, gstats] = fit(binned_f,binned_residuals,gamma_estimator, gamma_options);
    
    % <DEBUG>
    %figure
    %plot(binned_f, binned_residuals);
    %hold;
    %plot(gamma_fit);
    % </DEBUG>
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 4: Weigthed Least Square
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wtheta_weights = feval(gamma_fit,f);
    theta_options.StartPoint = NLLS_coeff;
    theta_options.Weights = wtheta_weights;

    [NLWLS_fit, wgof, wstats]= fit(f,T,theta_estimator, theta_options);
    if(sqrt(wgof.rmse)/norm(T)>0.001)
        warning('Second fit failed!');
    end
    NLWLS_coeff = coeffvalues(NLWLS_fit);
    % <DEBUG>
    %figure
    %plot(f,T)
    %hold;
    %plot(NLWLS_fit)
    % 
    %figure
    %plot(f,(stats.residuals));
    %hold;
    %plot(f,(wstats.residuals));
    % </DEBUG>
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5: Return value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = NLWLS_coeff(2);
    f_r = NLWLS_coeff(3);

    function [] = setupOptions()
        theta_estimator = fittype('theta_1/(1+theta_2^2*(x/theta_3-theta_3/x)^2)+theta_4');
        theta_options = fitoptions(theta_estimator);
        %theta_options.Robust = 'LAR';
        theta_options.MaxFunEvals = 10000;
        theta_options.MaxIter = 1000000;
        theta_options.DiffMinChange = min(T)*1e-5;
        theta_options.TolFun = min(T)*1e-10;
        theta_options.TolX = 1e-12; 
        theta_options.Lower = [min(T)*0.1 100 min(f) min(T)*0.1];
        theta_options.Upper = [1.5 100000 max(f) 1.5];
    end
end


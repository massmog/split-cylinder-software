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

function [ out ] = find_zeros(resonator, f_handle, n_I, sweep_var, varargin)
% Get correct arguments from function handle name (Fallback option)
handle_info = functions(f_handle);
handle_name = handle_info.function;

if(~isempty(findstr(handle_name,'Zmat')))
    varargin = {varargin{:} 'm12'};
end

if(~isempty(findstr(handle_name,'Jmat')))
    varargin = {varargin{:} 'janezic'};
end

if(~isempty(findstr(handle_name,'Cmat')))
    varargin = {varargin{:} 'm12' 'cal'};
end

%%% 
SINGLE = 0;
NODISP = 0;
CAL = 0;
nVarargs = length(varargin);
for k = 1:nVarargs
    switch varargin{k}
        case 'single'
            SINGLE = 1;
        case 'nodisp'
            NODISP = 1;
        case 'cal'
            CAL = 1;
    end
end

out = {};

TolX = 1e-100;
TolNull = 1e-10;
%options = optimoptions('fsolve','Display','Iter','TolX',TolX); % Option to display output
options = optimoptions('fsolve','TolX',TolX,'Display','off'); % Option to display output

if(numel(sweep_var)==2)
    Z = @(f_i) f_handle(f_i);
    det_f = @(f_i) arrayfun(@(x) real(det(Z(x)))+imag(det(Z(x))),f_i);
    I_max = max(sweep_var);
    I_min = min(sweep_var);   
    
    %% Search in subintervals
    delta = (I_max-I_min)/n_I;
    x_k = [I_min:delta:I_max];
    det_f_k = sign(det_f(x_k));

    for k = 1:n_I
       if(det_f_k(k)~=det_f_k(k+1))
            [x_0,fval,exit] = fzero(det_f,[x_k(k) x_k(k+1)],options);
            if(exit==1)
                %% Check Zero
                % Get null space of Z matrix
                [ Z_x_0, coefficients] = Z(x_0);
                %Using the null space command
                null_Z = null(Z_x_0);
                
                if isempty(null_Z) % Empty null of Z implies to solution
                    error_msg(sprintf('Invalid zero at\t %.3e - Zero deleted from output vector, no solution',x_0),varargin{:});
                elseif (sum(abs(null_Z)>1e-10)<=1) % Count valid elements, number of valid mode must be greater than 0
                    error_msg(sprintf('Invalid zero at\t %.3e - Zero deleted from output vector, no large elements',x_0),varargin{:});
                else
                    % Check norm of solution
                    nullity=norm(Z_x_0*null_Z,1);
                    % Mode identification and odd-mode check
                    A_n = null_Z(1:coefficients.N_u);
                    if(CAL)
                        B_n = null_Z(coefficients.N_u+1:coefficients.N_u);
                    else
                        B_n = null_Z(coefficients.N_u+1:coefficients.N_u+coefficients.N_s);
                    end
                    [A_max, A_max_n]=max(A_n);
                    
                    if(nullity>TolNull)
                        error_msg(sprintf('Invalid zero at \t %.3e - Null space check failed!',x_0),varargin{:});
                    elseif(abs(A_max)<0.3)
                        error_msg(sprintf('Invalid zero at \t %.3e - Mode is not dominant! (%e)',x_0,abs(A_max)),varargin{:});
                    else
                                       
                        %% Build identifier string
                        mode_name = sprintf('TE-0,%i',A_max_n);
                        
                        if(CAL==0)
                            [E_phi_max, H_z_u_loop, n_mode,E_phi_specimen]=plot_mode(resonator, coefficients, null_Z, mode_name, varargin{:});
                            if(abs(E_phi_specimen/E_phi_max)<1/4)
                                error_msg(sprintf('Invalid zero at\t %.3e - norm(Zx): %.1e - Substrate field strength is low!',x_0,nullity),varargin{:});
                                continue;
                            end
                        else
                            [E_phi_max, H_z_u_loop, n_mode]=plot_mode(resonator, coefficients, null_Z, mode_name, varargin{:});
                        end
                        
                        % Remove resonances with insufficient coupling
                        if(20*log10(abs(H_z_u_loop))<-10)
                            error_msg(sprintf('Invalid zero at\t %.3e - norm(Zx): %.1e - Insufficient coupling (<-10dB)!',x_0,nullity),varargin{:});
                            continue;
                        end
                        
                        error_msg(sprintf('Valid zero at\t %.3e - norm(Zx): %.1e',x_0,nullity),varargin{:})
                        
                        if(SINGLE)
                            out = x_0;
                            break;
                        end

                        if(~exist('out','var'))
                            out = {x_0;[A_max_n,n_mode];coefficients;null_Z};
                        else
                            out(:,size(out,2)+1) = {x_0;[A_max_n, n_mode];coefficients;null_Z};
                        end
                        
                    end                                  
                end
            end
        end
    end
else
    error('Error: Only one variable can be swept over a certain range. A maximum of two elements for start and end point of intervall are allowed.')
end

end


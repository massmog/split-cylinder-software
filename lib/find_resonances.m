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

function [f_r_Q] = find_resonances(f, S21,varargin)
    finder      = 0;
    plot_peaks  = 0;
    method      = '3dB';

    nVarargs = length(varargin);
    for k = 1:nVarargs
        switch varargin{k}
            case 'finder'
                finder = 1; 
            case 'plot_peaks'
                plot_peaks = 1;
            case 'estin'
                method  = 'estin';
            case 'coakley'
                method  = 'coakley';
            case 'circlefit'
                method  = 'circlefit';
        end
    end

    ma_offset   = 0;
    S21_finder  = S21;
    % Filter data for improved error resilience
    if(finder)
       % Moving average filter
       % Highest Q resonant frequency has lowest bandwidth = upper
       % bound for moving average

       Q_max    = 14000;
       f_max    = 20e9;
       d_f      = f_max/Q_max;

       n_points = length(f);
       f_range  = f(end)-f(1);

       % Safety factor 1/2 included!
       delta_n  = round(1/2*d_f/(f_range/n_points));

       windowSize =  delta_n;

       b = (1/windowSize)*ones(1,windowSize);
       S21_finder = filter(b,1,S21);

       % Moving average offset
       ma_offset = 10;
       S21_finder = S21_finder(ma_offset+1:end);
    end
    [S21_peaks,k_f] = findpeaks(abs(S21_finder),'MinPeakHeight',abs(max(S21))/15,'MinPeakProminence',abs((max(S21)-min(S21)))/15);
    k_f             = k_f + ma_offset;
    [~,k_f_min] = min(abs(f(k_f)-(f(end)+f(1))/2));

    if(plot_peaks)
        figure
        hLine = plot(f(ma_offset+1:end),20*log10(abs(S21_finder)),'Tag','Signal');
        hAxes = ancestor(hLine,'Axes');
        % turn on grid
        grid on;
        if numel(f)>1
          hAxes.XLim = hLine.XData([1 end]);
        end

        % use the color of the line
        color = get(hLine,'Color');
        hLine = line(hLine.XData(k_f-ma_offset),20*log10(abs(S21(k_f-ma_offset))),'Parent',hAxes, ...
             'Marker','o','LineStyle','none','Color',color,'tag','Peak');
    end

    k_range = zeros(length(k_f),2);
    for k = 1:length(k_f)
        if(k~=1)
            k_range(k,1)=round(mean([k_f(k-1) k_f(k)]));
        else
            k_range(k,1)=1;   
        end

        if(k~=length(k_f))
            k_range(k,2)=round(mean([k_f(k) k_f(k+1)]));
        else
            k_range(k,2)=length(f);
        end
    end
    % Resonance indicator
    res_ind = [0 mean_vec(diff(abs(S21)),round(0.1*length(S21)))];
    % If signal has too many maxima (noise!), do not record
    if(length(k_f)<10)
        f_r_Q = zeros(length(k_f),3);
        for k = 1:length(k_f)
            range = k_range(k,1):k_range(k,2);
            % Clip resonance curve to expected minimum Q 
            Q_min = 200;
            d_f = 5*f(k_f(k))/Q_min/2;
            %[~,n_max] = min(abs(f(range)-f(k_f(k))-d_f));
            %[~,n_min] = min(abs(f(range)-f(k_f(k))+d_f));
            % Clip resonance curve according to derivative
            d_k = round(0.1*d_f/(f(2)-f(1)));
            % Upper bound
            u_bnd = k_f(k)+d_k;
            if(u_bnd<=range(end))
                if(res_ind(u_bnd)<1)
                    b=u_bnd-1+find(res_ind(u_bnd:range(end))>0,1,'first');
                    if(~isempty(b))
                        range = range(1):b;
                    end
                end
            end
            
            % Lower bound
            l_bnd = k_f(k)-d_k;
            
            if(l_bnd>range(1))
                if(res_ind(l_bnd)<1)
                    b=range(1)-1+find(res_ind(range(1):l_bnd)>0,1,'last');
                    if(~isempty(b))
                        range = b:range(end);
                    end
                end
            end           

            % Clip range
            %range = range(n_min:n_max);
            
            switch method
                case '3dB'
                    [ f_r_Q(k,1), f_r_Q(k,2)]       = ThreeDBmethod( f(range), S21(range) );
                case 'estin'
                    [ f_r_Q(k,1), f_r_Q(k,2), ~]    = Estin( f(range), S21(range) );
                case 'coakley'
                    [ f_r_Q(k,1), f_r_Q(k,2)]       = Coakley( f(range), S21(range) );
                case 'circlefit'
                    if(k_f_min==k)
                        [ f_r_Q(k,1), f_r_Q(k,2)]       = CircleFit( f(range), S21(range), 'plot_fit');
                    else
                        [ f_r_Q(k,1), f_r_Q(k,2)]       = CircleFit( f(range), S21(range));
                    end
            end
            f_r_Q(k,3) = 20*log10(S21_peaks(k));
        end
    else
        f_r_Q = [];
    end
    
    function [avg_vec] = mean_vec(noise_vec,wSize)
        tabs = (1/wSize)*ones(1,wSize);
        avg_vec = filter(tabs,1,noise_vec);        
    end
end
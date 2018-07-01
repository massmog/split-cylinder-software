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

function [f_r_m, Q_m, a_l_m, sigma_m, cal_table] = measure_cal(f_max, cal_name,visaObj, varargin)
    varargin = {'m12' 'cal'};

    if(~exist('f_max','var'))
        f_max = 20e9;
    end
    
    if(~exist('cal_name','var'))        
        SAVE = 0;
    else
        SAVE = 1;
        % Define run number of this measurement
        run_no = -1;
        mm_folder = fullfile(pwd,'calibrations',sprintf('%s-%s',date,cal_name));
        run_file = fullfile(mm_folder,'run_id');
        
        if(~exist(mm_folder,'dir'))
            mkdir(mm_folder);
        end   
        
        if(exist(run_file,'file'))
            run_no = csvread(run_file);
        end
        
        run_no = run_no + 1;
        csvwrite(run_file,run_no);
    end


    %% Necessary includes
    addpath('./scsoft_m12')
    addpath('./lib');
    addpath('./qfactor');
    %% Constants
    resonator = constants();  

    %% Calculate Odd TE0n modes of resonator based on model
    zer = find_zeros(resonator, @(x) Cmat(resonator, x),1000,[10e9 f_max],varargin{:})
    
    %% Measure all resonances
    % SETUP 
    Q_min = 500;
    
    n_zer = size(zer,2);
    f_r_m = zeros(1,n_zer);
    Q_m =  zeros(1,n_zer);
    
    a_l_m = zeros(1,n_zer);
    sigma_m = zeros(1,n_zer);
    
    cal_table = zeros(4,n_zer);
    
    %% SERIAL COMMUNICATION
    ser = setupSerial();    
    %% GPIB Communication
    %setupGPIB();    
    
    for l = 1:n_zer
        f_r = zer{1,l};
        d_f = f_r/Q_min/2;
        
        %%<RELEASE>
        % Find exact resonance frequency around estimate
        % First test measurement
        measurement_range = [f_r-d_f,f_r+d_f];
        [f,S21] = visa_S21(measurement_range,10000,1,visaObj);
        % Sort resonances by Q
        [f_0,~,S] =  get_resonance(f,S21,'finder');  
        
        if(~isempty(f_0))       
          
            k = 0;
            while((S<-65)||(S>-55))
                n_steps = 5;

                % FEED ADJUSTMENTS
                if(S>-55)
                    up_feed(n_steps);
                end
                if(S<-65)
                    down_feed(n_steps);
                end

                % TIMEOUT
                k = k+1;
                if(k>round(330/n_steps))
                    break;
                end

                [f,S21] = visa_S21(measurement_range,20001,1,visaObj,'no_setup');
                [~,~,S] =  get_resonance(f,S21,'finder');
            end
        
      
            [f,S21] = visa_S21(measurement_range,20001,25,visaObj);
            [f_r_m(l),Q_m(l),S] =  get_resonance(f,S21,'plot_peaks','circlefit');
            
        else
            continue;
        end
        %%</RELEASE>
        
        %%<DEBUG>
        %f_r_m(l) = zer{1,l};
        %Q_m(l) = 1000;
        %%</DEBUG>
    
        %%% Error handling!!! No zeros found....
        
        
        %% All resonances measured, Calibration procedure
        cal_zer = find_zeros(resonator, @(x) Cmat(resonator, f_r_m(l),'a_l',x),100,[resonator.a_l-1e-4 resonator.a_l+1e-4],varargin{:});
        a_l_m(l) = cal_zer{1,1};
        resonator.a_l = a_l_m(l);
        coefficients = cal_zer{3,1};
        null_Z = cal_zer{4,1};
        mn_mode = cal_zer{2,1}';
        sigma_m(l) = Closs( resonator, coefficients, null_Z, 'sigma',Q_m(l));
        
        cal_table(:,l)=[mn_mode,a_l_m(l),sigma_m(l)];
        
        % Save measured resonance curves
        if(SAVE)
                measurement_file = fullfile(mm_folder,sprintf('mode-%i-%i.mat',run_no,l)); 
                save(measurement_file,'f','S21')
        end
        
    end
    % Save measurement results
    if(SAVE)
        result_file = fullfile(mm_folder,sprintf('result-%i.mat',run_no));
        save(result_file,'f_r_m','Q_m','a_l_m','sigma_m');
    end
    
    % Close serial connection
    destroySerial();
    % Close VISA object
    %destroyGPIB();
    
    function [] = up_feed(N_loop)
        for g = 1:N_loop
            fprintf(ser,'u');
            read_str = '';
            while(~isequal(read_str,'rdy'))
                %NOP 
                read_str = fscanf(ser);
                read_str = read_str(isletter(read_str));
                pause(1e-6);
            end
        end
    end

    function [] = down_feed(N_loop)
        for g = 1:N_loop
            fprintf(ser,'d');
            read_str = '';
            while(~isequal(read_str,'rdy'))
                %NOP 
                read_str = fscanf(ser);
                read_str = read_str(isletter(read_str));
                pause(1e-6);
            end 
        end
    end

    function [f,Q,S] = get_resonance(f, S21, varargin)
        try
            [f_Q] =  sortrows(find_resonances(f,S21,varargin{:}),2);
        catch
            warning('Resonance could not be found!');
            return;
        end
        
        if(~isempty(f_Q))
            f = f_Q(end,1);
            Q = f_Q(end,2);
            S = f_Q(end,3);
        else
            f = [];
            Q = [];
            S = [];
        end
    end

    function [ser] = setupSerial()
            % Clean up
            oldObjects=instrfind;
            if ~isempty(oldObjects)
                for k = 1:length(oldObjects)
                    if(isequal(oldObjects(k).type,'serial'))
                        delete(oldObjects(k));
                        clear oldObjects(k);
                    end
                end
            end    
            % Start serial connection    
            ser = serial('com7');
            fopen(ser);
    end
    
    function [] = destroySerial()
            fclose(ser);
            delete(ser);
            clear ser;
    end

    function [] = setupGPIB()
            % Clean up
         oldObjects=instrfind;
        if ~isempty(oldObjects)
            for o = 1:length(oldObjects)
                if(isequal(oldObjects(o).type,'visa-gpib'))
                    delete(oldObjects(o));
                    clear oldObjects(o);
                end
            end
        end

        % Start GPIB connection through VISA
        resourceStruct = instrhwinfo('visa','agilent'); %#ok<NASGU>
        % A generalized way to create a visa interface object:
        % eval(['visaObj = ',res.ObjectConstructorName{1}]);
        visaObj = visa('agilent','GPIB0::16::INSTR');
        
        numOfPoints = 20001;
        %% Configure interface object
        % Set a sufficiently large input buffer size to store the S-Parameter data
        buffer_size = numOfPoints*72;
        set(visaObj, 'InputBufferSize', buffer_size);
        % Set large timeout in the event of long s-parameter measurement
        set(visaObj, 'Timeout', 600);
        
        fopen(visaObj);
    end
    
    function [] = destroyGPIB()
            delete(visaObj);
            clear visaObj;
    end
end
 

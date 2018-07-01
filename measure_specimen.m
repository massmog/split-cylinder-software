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

function [f_r,Q, e_r, tand] = measure_specimen(d, f_max, measurement_name,visaObj)
    if(~exist('f_max','var'))
        f_max = 20e9;
    end
    
    if(~exist('measurement_name','var'))        
        SAVE = 0;
    else
        SAVE = 1;
        % Define run number of this measurement
        run_no = -1;
        mm_folder = fullfile(pwd,'measurements',sprintf('%s-%s',date,measurement_name));
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

    %% SERIAL COMMUNICATION
    ser = setupSerial();   
    
    %% GPIB Communication
    %setupGPIB();
    
    % Find fundamental mode
    % TE-11 resonance frequency for the empty resonator f_r_TE11 =
    % 5.4950e+09 - Maximum frequency limited to 5.6GHz
    [f,S21] = visa_S21([4e9,5.6e9],20001,5,visaObj);
    % Estimate permittivity from fundamental TE111 mode
    f_r_Q_TE111 = find_resonances(f,abs(S21),'plot_peaks','finder');
    f_r_TE111 = f_r_Q_TE111(1,1);
    e_r_est = TE_estimate(1,1,resonator,f_r_TE111,d);
      
    % Calculate permittivity and loss factor from detail measurement 
    e_r_last = e_r_est; 
    f_r_last = 1e9;
    d_f = 0;
    
    f_r = [];
    Q = [];
    e_r = [];
    tand = [];
    
    for m = 1:100
        next_higher_mode = find_zeros(resonator,@(x) Zmat(resonator,e_r_est,x,d),1000,[f_r_last+d_f/5 f_r_last+10e9],'nodisp','single');
        
        f_r_est = next_higher_mode;
        if(f_r_est>f_max)
            break;
        end
        % TO-DO: Make this dependent von different parameters (thickness,
        % e_r...) % Include coupling handling
        Q_min = 200;
        % TO-DO: Self-adjusting width, zoom function
        d_f = 5*f_r_est/Q_min/2;
        measurement_range = [f_r_est-d_f,f_r_est+d_f];
        
        [f,S21] = visa_S21(measurement_range,5000,10,visaObj);    
        f_Q_S = find_resonances(f,abs(S21),'finder'); 
        %% TO-DO Plausibility check, if multiple resonances use the one with the most promising values
        % Choose resonance with smallest distance from measured frequency
        
        %% TO-DO INCLUDE LOW-POWER COMPENSATION
        if(~isempty(f_Q_S))
            [~,r_min] = min(abs(f_Q_S(:,1)-f_r_est));
            S = f_Q_S(r_min,3);
            
            k = 0;
            S_u_lim = -55;
            S_l_lim = -65;
            while((S<S_l_lim)||(S>S_u_lim))
                n_steps = 5;

                % FEED ADJUSTMENTS
                if(S>S_u_lim)
                    up_feed(n_steps);
                end
                if(S<S_l_lim)
                    down_feed(n_steps);
                end

                % TIMEOUT
                k = k+1;
                if(k>round(330/n_steps))
                    break;
                end

                [f,S21] = visa_S21(measurement_range,5000,5,visaObj,'no_setup');
                f_Q_S = find_resonances(f,abs(S21),'finder'); 
                if(~isempty(f_Q_S))
                    [~,r_min] = min(abs(f_Q_S(:,1)-f_r_est));
                    S = f_Q_S(r_min,3);
                else
                    break;
                end
            end
        
      
            [f,S21] = visa_S21(measurement_range,20001,25,visaObj);
            f_Q_S = find_resonances(f,S21,'circlefit');            
            
            if(~isempty(f_Q_S))
                [~,r_min] = min(abs(f_Q_S(:,1)-f_r_est));
                f_r = [f_r,f_Q_S(r_min,1)];
                Q = [Q,f_Q_S(r_min,2)];

                % Compute permittivity
                zer = find_zeros(resonator,@(x) Zmat(resonator,x,f_Q_S(r_min,1),d),1000,[e_r_last*0.8 1.3*e_r_last]);
                % If multiple solutions are found the solution closest to the
                % previous e_r value is taken as the solution
                [~,e_min]=min(abs([zer{1,:}]-e_r_last));

                e_r = [e_r,zer{1,e_min}]
                e_r_last = zer{1,e_min};

                % Compute tangent delta
                coeff = zer{3,e_min};
                null_Z = zer{4,e_min};

                tand=[tand,Zloss( resonator, coeff, null_Z,'tand',f_Q_S(r_min,2))]
            end
        end
        
        % Save measured resonance curves
        if(SAVE)
            measurement_file = fullfile(mm_folder,sprintf('mode-%i-%i.mat',run_no,m)); 
            save(measurement_file,'f','S21')
        end
        f_r_last = f_r_est;
    end
    
    % Save measurement results
    if(SAVE)
        result_file = fullfile(mm_folder,sprintf('result-%i.mat',run_no));
        save(result_file,'f_r','Q','e_r','tand','d','f_r_Q_TE111');
    end
    
    % Close serial connection
    destroySerial();
    % Close VISA object
    %destroyGPIB();
    
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
        set(visaObj, 'Timeout', 60);
        
        fopen(visaObj);
    end

    function [] = destroyGPIB()
        fclose(visaObj);
        delete(visaObj);
        clear visaObj;
    end
%%% END SUBROUTINES %%%
    
end

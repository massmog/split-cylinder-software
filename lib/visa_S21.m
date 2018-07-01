Copyright (c) 2016, The MathWorks, Inc. 
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are 
met:

* Redistributions of source code must retain the above copyright 
notice, this list of conditions and the following disclaimer. 
* Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimer in 
the documentation and/or other materials provided with the distribution. 
* In all cases, the software is, and all modifications and derivatives 
of the software shall be, licensed to you solely for use in conjunction 
with MathWorks products and service offerings.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

function [frequencies S21] = visa_S21(measurement_range,numOfPoints,num_avg, visaObj,varargin)
    no_setup = 0;
    view_only = 0;
    plot_SNP = 0;
    
    nVarargs = length(varargin);
    for k = 1:nVarargs
        switch varargin{k}
            case 'no_setup'
                no_setup = 1; 
            case 'view_only'
                view_only = 1;
            case 'plot'
                plot_SNP = 1;
        end
    end
    
    f_start = measurement_range(1);
    f_stop = measurement_range(2);
    
    %% Baseband Equivalent Modeling of a Transmission Line

    % In this example, you model an RF transmission line stimulated by a pulse
    % and plot the baseband-equivalent model that the blockset uses to simulate
    % the transmission line in the time domain. This example helps you 
    % understand how to best apply the baseband-equivalent modeling paradigm of
    % performing time-domain simulation using a limited band of frequency data.
    % This demo uses SCPI commands to interface with a PNA Series Network
    % Analyzer (from Agilent). It is assumed that:
    %
    % 1. Calibration of the PNA has been completed by the user.
    % 2. The Agilent IO libraries have been installed. 
    % 
    % Author: Siddhartha Shankar

    if(~exist('numOfPoints'))
        numOfPoints = 201;
    end   
    
    %% Demo Requirements
    %
    % This demonstration requires the following products:
    %
    %    * MATLAB       
    %    * Instrument Control Toolbox
    %
    % The Following products are required to utilize the optional sections of
    % this demo:
    %
    %   * Simulink
    %   * RF Toolbox
    %   * RF Blockset
    %% Additional Notes
    %
    % * This M-File was last used and tested with MATLAB 7.6 (R2008a)
    % * Consider using the MATLAB Report Generator to create a report from this
    %   M-File

    %%% Configure interface object
    %% Set a sufficiently large input buffer size to store the S-Parameter data
    %buffer_size = numOfPoints*72;
    %set(visaObj, 'InputBufferSize', buffer_size);
    %% Set large timeout in the event of long s-parameter measurement
    %set(visaObj, 'Timeout', 60);
    
    %fopen(visaObj);

    %% Configure parameter to be acquired and initiate 
    if(~no_setup)
        localConfigurePNA('MA');
        setupMeasurement();
        setupTrace();
    end
    if(~view_only)
        data = localFetchData(visaObj,'S21');
    end
    clrdevice(visaObj);
    % Format of returned data is as follows:
    % Row 1 - Frequency: f1 f2 f3 ... fn
    % Row 2 - Smn Magnitude: MX1 MX2 MX3 ... MXn (m = 1 to 2)
    % Row 3 - Smn Angle: AX1 AX2 AX3 ... AXn (n = 1 to 2)

    % Read frequency data back from returned data
    if(~view_only)
        frequencies=data(1,:);
        magS21=data(4,:);
        degS21=data(5,:);
        S21 = magS21.*exp(i*deg2rad(degS21));
        if(plot_SNP)
            figure
            plot(frequencies,20*log10(abs(S21)));
            title('S21 mag');
            xlabel('f');
            ylabel('dB');
        end
    end
    %% Close connection, delete visa object
    %fclose(visaObj);

    function localWaitForSystemReady(visaObj)
        opcStatus = 0;
        while(~opcStatus)
            opcStatus = query(visaObj, '*OPC?','%s\n','%d'); 
        end
    end
    function localWaitForAveragingFinished(visaObj)
        avgStatus = 0;
        while(~avgStatus)
            avgStatus = (query(visaObj, 'STAT:OPER:AVER1:COND?','%s\n','%i')==2); 
        end
    end    
    function localConfigurePNA(fileFormat)
        %% Preset system
        % SYSTem:PRESet
        % OPC? = All Operations Complete? +1 for Yes
        fprintf(visaObj,'SYST:PRES');
        localWaitForSystemReady(visaObj);

        %% Set S2P File Format. 
        fprintf(visaObj, sprintf('MMEM:STOR:TRAC:FORM:SNP %s',fileFormat)); 
        % MA - Linear Magnitude / degrees
        % DB - Log Magnitude / degrees
        % RI - Real / Imaginary
        % AUTO - data is output in currently selected trace form

        %% Set byte order to swapped (little-endian) format
        % FORMat:BORDer <char>
        fprintf(visaObj, 'FORM:BORD SWAP');
        % NORMal - Use when your controller is anything other than an IBM compatible computers
        % SWAPped - for IBM compatible computers
        
        %% 
        %fprintf(visaObj, 'DISP:WIND1:STATE ON');
        %% Set data type to real 64 bit binary block
        % FORMat[:DATA] <char>, 64 for more significant digits and precision
        fprintf(visaObj, 'FORM REAL,64');
        % REAL,32 - (default value for REAL) Best for transferring large amounts of measurement data.
        % REAL,64 - Slower but has more significant digits than REAL,32. Use REAL,64 if you have a computer that doesn't support REAL,32.
        % ASCii,0 - The easiest to implement, but very slow. Use if small amounts of data to transfer.
    end
    function data = localFetchData(visaObj,sParameter)
        % Set up the trace corresponding to PARAMETER on the PNA and return DATA,
        % a matrix of 2-port S-Parameters in S2P format with specified PRECISION.
        % COUNT is the number of values read and MESSAGE tells us if the read
        % operation was unsuccessful for some reason.
        fprintf(visaObj,sprintf('CALC:PAR:MOD %s',sParameter));
        localWaitForSystemReady(visaObj);
        localWaitForAveragingFinished(visaObj);
        fprintf(visaObj, 'CALC:DATA:SNP? 2');
        %localWaitForSystemReady(visaObj);
        % Read the data back using binblock format
        [rawData] = binblockread(visaObj, 'double');
        data = reshape(rawData, [(length(rawData)/9),9]);
        data = data';
    end
    function [] = setupMeasurement()
        % Visual Confirmation of S21 in DB format. 
        % Select S21
        % Delete all measurements on the PNA
        fprintf(visaObj,'CALC:PAR:DEL:ALL');
        % Create S21 measurement "MyMeas"
        fprintf(visaObj,'CALC:PAR:DEF:EXT "MyMeas",S21');
        % Add trace to display
        fprintf(visaObj,'DISP:WIND1:TRAC1:FEED "MyMeas"');
        localWaitForSystemReady(visaObj);
        % Set output power
        fprintf(visaObj,'SOUR:POW1 0');
        localWaitForSystemReady(visaObj);
        
        if(exist('num_avg','var'))
            fprintf(visaObj,'SENS1:AVER:STAT ON');
            fprintf(visaObj,sprintf('SENS1:AVER:COUN %i',num_avg));
        end
        fprintf(visaObj,sprintf('SENS:FREQ:STAR %i',f_start));
        localWaitForSystemReady(visaObj);
        fprintf(visaObj,sprintf('SENS:FREQ:STOP %i',f_stop));
        localWaitForSystemReady(visaObj);
        fprintf(visaObj,sprintf('SENS:SWE:POIN %i',numOfPoints));
        localWaitForSystemReady(visaObj);
    end
    function [] = setupTrace()
        fprintf(visaObj,'DISP:WIND:TRAC:Y:PDIV 10');
        localWaitForSystemReady(visaObj);
        fprintf(visaObj,'DISP:WIND:TRAC:Y:RLEV -80');
        localWaitForSystemReady(visaObj);
        fprintf(visaObj,'DISP:WIND:TRAC:Y:RPOS 4'); 
        localWaitForSystemReady(visaObj);
        fprintf(visaObj, 'CALC:PAR:SEL "MyMeas"'); 
        localWaitForSystemReady(visaObj);
    end
end
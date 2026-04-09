classdef ErrorMetric < handle
    properties
        name          % e.g., 'SIR', 'ENKF'
        DxArray       % vector of Dx values
        NR            % number of repetitions
        MSEf          % numel(DxArray)xNR matrix
        MSEp          % same
        runtime       %% Execution time in seconds (Dx x NR)
    end
    
    methods
        function obj = ErrorMetric(name, DxArray, NR)
            obj.name = name;
            obj.DxArray = DxArray;
            obj.NR = NR;
            nDx = numel(DxArray);
            obj.MSEf = nan(nDx, NR);
            obj.MSEp = nan(nDx, NR);
            obj.runtime = nan(nDx, NR);
        end
        
        function record(obj, iDX, nr, MSEf, MSEp,time_sec)
            obj.MSEf(iDX, nr) = MSEf;
            obj.MSEp(iDX, nr) = MSEp;
            obj.runtime(iDX, nr) = time_sec;
        end
        
        function disp(obj)
            disp(table( ...
                obj.DxArray', ...
                mean(obj.MSEf,2), ...
                mean(obj.MSEp,2), ...
                mean(obj.runtime,2), ...
                'VariableNames', {'Dx','MSEf_mean','MSEp_mean','Runtime_sec'} ...
            ));
        end
    end
end

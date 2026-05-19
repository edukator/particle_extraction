classdef ErrorMetric < handle
    properties
        name          % e.g., 'SIR', 'ENKF'
        DxArray       % vector of Dx values
        MSEf          % numel(DxArray)x1 vector
        MSEp          % numel(DxArray)x1 vector
        runtime       % Execution time in seconds (Dx x 1)
    end
    
    methods
        function obj = ErrorMetric(name, DxArray)
            obj.name = name;
            obj.DxArray = DxArray;
            nDx = numel(DxArray);
            obj.MSEf = nan(nDx, 1);
            obj.MSEp = nan(nDx, 1);
            obj.runtime = nan(nDx, 1);
        end
        
        function record(obj, iDX, MSEf, MSEp,time_sec)
            obj.MSEf(iDX) = MSEf;
            obj.MSEp(iDX) = MSEp;
            obj.runtime(iDX) = time_sec;
        end
        
        function disp(obj)
            disp(table( ...
                obj.DxArray', ...
                obj.MSEf, ...
                obj.MSEp, ...
                obj.runtime, ...
                'VariableNames', {'Dx','MSEf_mean','MSEp_mean','Runtime_sec'} ...
            ));
        end
    end
end

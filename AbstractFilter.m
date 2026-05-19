classdef (Abstract) AbstractFilter < handle
    properties
        params  % struct of shared simulation data
        metric   % ErrorMetric object
    end
   
    methods
        function obj = AbstractFilter(params, name, DxArray)
            obj.params = params;
            obj.metric = ErrorMetric(name, DxArray);
        end
        function run(obj, iDX)
            tStart = tic;
            [Xf, Xp] = obj.runFilter();
            elapsed = toc(tStart);

            p = obj.params;
            if any(~isfinite(Xf(:)))
                MSEf = NaN;
                MSEp = NaN;
                obj.metric.record(iDX, MSEf, MSEp, elapsed);
                return;
            end
            MSEf = mean(sum((Xf - p.truth).^2)) / p.Pd_f;
            MSEp = mean(sum((Xp - p.truth_full).^2)) / p.Pd_p;

            obj.metric.record(iDX, MSEf, MSEp, elapsed);
        end
    end
     methods (Abstract, Access = protected)
        [Xf, Xp] = runFilter(obj)
    end
end

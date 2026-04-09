classdef SIRFilter < AbstractFilter
    methods
        function obj = SIRFilter(params, DxArray, NR)
            obj@AbstractFilter(params, 'SIR', DxArray, NR);
        end
    end
     methods (Access = protected)

        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp] = sir(p.F, p.sx, p.sz, ...
                          p.he, p.NTe, p.n_obs, ...
                          p.ze_sparse, p.H, p.X0, p.ness_thr);
        end
    
    end
end

classdef APFFilter < AbstractFilter
    methods
        function obj = APFFilter(params, DxArray, NR)
            obj@AbstractFilter(params, 'APF', DxArray, NR);
        end
    end
    methods (Access = protected)
        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp] = apf_new(p.F, p.sx, p.sz, ...
                               p.he, p.NTe, p.n_obs, ...
                               p.ze_sparse, p.H, p.X0);
        end
    end
end

classdef ENKFFilter < AbstractFilter
    methods
        function obj = ENKFFilter(params, DxArray, NR)
            obj@AbstractFilter(params, 'EnKF', DxArray, NR);
        end
    end
    methods (Access = protected)

        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp] = enkfH_e(p.F, p.sx, p.sz, ...
                              p.he, p.NTe, p.n_obs, ...
                              p.ze_sparse, p.H, p.X0, ...
                              p.Dz, p.fixed_observed_components);
        end
    end
end

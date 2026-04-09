classdef ENKF_Barrier_Filter < AbstractFilter
    methods
        function obj = ENKF_Barrier_Filter(params, DxArray, NR)
            obj@AbstractFilter(params, 'Barrier-EnKF', DxArray, NR);
        end
    end
    methods (Access = protected)

        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp] = enkfH_e_Barrier(p.F, p.sx, p.sz, ...
                              p.he, p.NTe, p.n_obs, ...
                              p.ze_sparse, p.H, p.X0, ...
                              p.Dz, p.fixed_observed_components,p.barrier_params);
        end
    end
end

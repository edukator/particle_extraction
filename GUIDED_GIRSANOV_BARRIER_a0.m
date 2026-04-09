classdef GUIDED_GIRSANOV_BARRIER_a0 < AbstractFilter
    methods
        function obj = GUIDED_GIRSANOV_BARRIER_a0(params, DxArray, NR)
            obj@AbstractFilter(params, 'Guided_Girsanov_Barrier_a0', DxArray, NR);
        end
    end
    methods (Access = protected)

        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp]= sir_guided_girsanov_barrier_a0(p.F,p.sx,p.sz,p.he,p.NTe,p.n_obs, ...
                p.ze_sparse,p.H,p.X0,p.ness_thr, p.barrier_params);
            
        end
    end
end

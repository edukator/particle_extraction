classdef BarrierFilter < AbstractFilter
    methods
        function obj = BarrierFilter(params, DxArray, NR)
            obj@AbstractFilter(params, 'Barrier-SIR', DxArray, NR);
        end
    end
    methods (Access = protected)
        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp,dummy] = sir_barrier(p.F, p.sx, p.sz, ...
                                      p.he, p.NTe, p.n_obs, ...
                                      p.ze_sparse, p.H, p.X0, ...
                                      p.ness_thr,...
                                       p.barrier_params);
        end
    end
end

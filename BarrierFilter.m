classdef BarrierFilter < AbstractFilter
    methods
        function obj = BarrierFilter(params, DxArray)
            obj@AbstractFilter(params, 'Barrier-SIR', DxArray);
        end
    end
    methods (Access = protected)
        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            save_obs_indices = [];
            if isfield(p, 'save_obs_indices')
                save_obs_indices = p.save_obs_indices;
            end
            snapshot_cfg = struct();
            if isfield(p, 'snapshot_dir')
                snapshot_cfg.dir = p.snapshot_dir;
                snapshot_cfg.filter_name = 'Barrier_SIR';
                snapshot_cfg.Dx = p.Dx;
            end
            [Xf, Xp] = sir_barrier(p.F, p.sx, p.sz, ...
                                      p.he, p.NTe, p.n_obs, ...
                                      p.ze_sparse, p.H, p.X0, ...
                                      p.ness_thr,...
                                       p.barrier_params, save_obs_indices, snapshot_cfg);
        end
    end
end

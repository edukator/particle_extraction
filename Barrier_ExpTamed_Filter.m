classdef Barrier_ExpTamed_Filter < AbstractFilter
    methods
        function obj = Barrier_ExpTamed_Filter(params, DxArray)
            obj@AbstractFilter(params, 'Barrier-ExpTamed-SIR', DxArray);
        end
    end
    methods (Access = protected)
        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            save_step_indices = [];
            if isfield(p, 'save_step_indices')
                save_step_indices = p.save_step_indices;
            elseif isfield(p, 'save_obs_indices')
                save_step_indices = p.save_obs_indices;
            end
            snapshot_cfg = struct();
            if isfield(p, 'snapshot_dir')
                snapshot_cfg.dir = p.snapshot_dir;
                snapshot_cfg.filter_name = 'Barrier-ExpTamed-SIR';
                snapshot_cfg.Dx = p.Dx;
                snapshot_cfg.fixed_observed_components = p.fixed_observed_components;
                snapshot_cfg.non_observed_components = setdiff(1:p.Dx, p.fixed_observed_components);
                snapshot_cfg.truth_full = p.truth_full;
                snapshot_cfg.ze_sparse = p.ze_sparse;
                snapshot_cfg.n_obs = p.n_obs;
                snapshot_cfg.H = p.H;
            end
            [Xf, Xp] = sir_barrier_tamed(p.F, p.sx, p.sz, ...
                                      p.he, p.NTe, p.n_obs, ...
                                      p.ze_sparse, p.H, p.X0, ...
                                      p.ness_thr,...
                                       p.barrier_params, save_step_indices, snapshot_cfg,'ExpTamed');
        end
    end
end

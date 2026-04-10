classdef ENKF_Barrier_Filter < AbstractFilter
    methods
        function obj = ENKF_Barrier_Filter(params, DxArray)
            obj@AbstractFilter(params, 'Barrier-EnKF', DxArray);
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
                snapshot_cfg.filter_name = 'Barrier_EnKF';
                snapshot_cfg.Dx = p.Dx;
            end
            [Xf, Xp] = enkfH_e_Barrier(p.F, p.sx, p.sz, ...
                              p.he, p.NTe, p.n_obs, ...
                              p.ze_sparse, p.H, p.X0, ...
                              p.Dz, p.fixed_observed_components,p.barrier_params, save_obs_indices, snapshot_cfg);
        end
    end
end

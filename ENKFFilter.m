classdef ENKFFilter < AbstractFilter
    methods
        function obj = ENKFFilter(params, DxArray)
            obj@AbstractFilter(params, 'EnKF', DxArray);
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
                snapshot_cfg.filter_name = 'EnKF';
                snapshot_cfg.Dx = p.Dx;
            end
            [Xf, Xp] = enkfH_e(p.F, p.sx, p.sz, ...
                              p.he, p.NTe, p.n_obs, ...
                              p.ze_sparse, p.H, p.X0, ...
                              p.Dz, p.fixed_observed_components, save_step_indices, snapshot_cfg);
        end
    end
end

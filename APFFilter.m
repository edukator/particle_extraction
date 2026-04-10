classdef APFFilter < AbstractFilter
    methods
        function obj = APFFilter(params, DxArray)
            obj@AbstractFilter(params, 'APF', DxArray);
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
                snapshot_cfg.filter_name = 'APF';
                snapshot_cfg.Dx = p.Dx;
            end
            [Xf, Xp] = apf_new(p.F, p.sx, p.sz, ...
                               p.he, p.NTe, p.n_obs, ...
                               p.ze_sparse, p.H, p.X0, save_obs_indices, snapshot_cfg);
        end
    end
end

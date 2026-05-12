classdef SIRFilter < AbstractFilter
    methods
        function obj = SIRFilter(params, DxArray)
            obj@AbstractFilter(params, 'SIR', DxArray);
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
                snapshot_cfg.filter_name = 'SIR';
                snapshot_cfg.Dx = p.Dx;
            end
            [Xf, Xp] = sir(p.F, p.sx, p.sz, ...
                          p.he, p.NTe, p.n_obs, ...
                          p.ze_sparse, p.H, p.X0, p.ness_thr, save_step_indices, snapshot_cfg);
        end
    
    end
end

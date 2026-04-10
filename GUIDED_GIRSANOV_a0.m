classdef GUIDED_GIRSANOV_a0 < AbstractFilter
    methods
        function obj = GUIDED_GIRSANOV_a0(params, DxArray)
            obj@AbstractFilter(params, 'Guided_Girsanov_a0', DxArray);
        end
    end
    methods (Access = protected)

        function [Xf, Xp] = runFilter(obj)
            p = obj.params;
            [Xf, Xp]= sir_guided_girsanov_a0(p.F,p.sx,p.sz,p.he,p.NTe,p.n_obs, ...
                p.ze_sparse,p.H,p.X0,p.ness_thr);
            
        end
    end
end

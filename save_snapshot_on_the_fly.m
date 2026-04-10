function save_snapshot_on_the_fly(snapshot_cfg, save_obs_indices, obs_idx, particles, weights)
% Save particle/weight snapshots to disk during filter execution.
%
% snapshot_cfg fields (optional):
%   - dir: output directory
%   - filter_name: string used in output filename
%   - Dx: state dimension (metadata)

    if nargin < 5
        return;
    end
    if isempty(save_obs_indices) || ~ismember(obs_idx, save_obs_indices)
        return;
    end
    if isempty(snapshot_cfg) || ~isstruct(snapshot_cfg)
        return;
    end
    if ~isfield(snapshot_cfg, 'dir') || isempty(snapshot_cfg.dir)
        return;
    end

    if ~exist(snapshot_cfg.dir, 'dir')
        mkdir(snapshot_cfg.dir);
    end

    filter_name = 'filter';
    if isfield(snapshot_cfg, 'filter_name') && ~isempty(snapshot_cfg.filter_name)
        filter_name = snapshot_cfg.filter_name;
    end

    snapshot = struct();
    snapshot.obs_idx = obs_idx;
    snapshot.particles = particles;
    snapshot.weights = weights;
    if isfield(snapshot_cfg, 'Dx')
        snapshot.Dx = snapshot_cfg.Dx;
    end

    filename = sprintf('%s_obs_%04d.mat', filter_name, obs_idx);
    save(fullfile(snapshot_cfg.dir, filename), 'snapshot');
end

function save_snapshot_on_the_fly(snapshot_cfg, save_step_indices, current_step, particles, weights)
    if nargin < 5
        return;
    end
    % Check if the current absolute step is in the list of steps to save
    if isempty(save_step_indices) || ~ismember(current_step, save_step_indices)
        return;
    end
    if isempty(snapshot_cfg) || ~isstruct(snapshot_cfg) || ~isfield(snapshot_cfg, 'dir') || isempty(snapshot_cfg.dir)
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
    snapshot.step_idx = current_step; % Store step instead of obs_idx
    snapshot.particles = particles;
    snapshot.weights = weights;
    if isfield(snapshot_cfg, 'Dx')
        snapshot.Dx = snapshot_cfg.Dx;
    end

    % Update filename to reflect absolute steps (e.g., SIR_step_00150.mat)
    filename = sprintf('%s_step_%05d.mat', filter_name, current_step);
    save(fullfile(snapshot_cfg.dir, filename), 'snapshot');
end
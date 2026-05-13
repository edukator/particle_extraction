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
    if isfield(snapshot_cfg, 'fixed_observed_components')
        snapshot.fixed_observed_components = snapshot_cfg.fixed_observed_components(:)';
    else
        snapshot.fixed_observed_components = [];
    end
    if isfield(snapshot_cfg, 'non_observed_components')
        snapshot.non_observed_components = snapshot_cfg.non_observed_components(:)';
    else
        snapshot.non_observed_components = [];
    end

    if isfield(snapshot_cfg, 'truth_full') && ~isempty(snapshot_cfg.truth_full)
        truth_idx = current_step + 1; % truth_full is stored from t=0
        if truth_idx >= 1 && truth_idx <= size(snapshot_cfg.truth_full, 2)
            snapshot.true_state = snapshot_cfg.truth_full(:, truth_idx);
        else
            snapshot.true_state = [];
        end
    else
        snapshot.true_state = [];
    end

    if isfield(snapshot_cfg, 'ze_sparse')
        snapshot.ze_sparse = snapshot_cfg.ze_sparse;
    end
    if isfield(snapshot_cfg, 'n_obs')
        snapshot.n_obs = snapshot_cfg.n_obs;
    end
    if isfield(snapshot_cfg, 'H')
        snapshot.H = snapshot_cfg.H;
    end

    if ~isempty(snapshot.fixed_observed_components)
        snapshot.observed_particles = particles(snapshot.fixed_observed_components, :);
    else
        snapshot.observed_particles = [];
    end
    if ~isempty(snapshot.non_observed_components)
        snapshot.non_observed_particles = particles(snapshot.non_observed_components, :);
    else
        snapshot.non_observed_particles = [];
    end

    % Update filename to reflect absolute steps (e.g., SIR_step_00150.mat)
    filename = sprintf('%s_step_%05d.mat', filter_name, current_step);
    save(fullfile(snapshot_cfg.dir, filename), 'snapshot');
end

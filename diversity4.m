% =========================================================================
% Compare Barrier-SIR and SIR particle densities in a paper-sized figure
% =========================================================================
barrier_file = 'sim1/snapshots/Dx_1000/Barrier_SIR_step_10000.mat';
sir_file = 'sim1/snapshots/Dx_1000/SIR_step_10000.mat';

barrier_snapshot = load_snapshot(barrier_file);
sir_snapshot = load_snapshot(sir_file);

if barrier_snapshot.Dx ~= sir_snapshot.Dx
    error('The snapshots must have the same state dimension.');
end
if barrier_snapshot.step_idx ~= sir_snapshot.step_idx
    error('The snapshots must have the same step index.');
end

Dx = barrier_snapshot.Dx;
step_idx = barrier_snapshot.step_idx;
he = 1e-3;
time_val = step_idx * he;

[barrier_observed, barrier_hidden] = get_components(barrier_snapshot);
[sir_observed, sir_hidden] = get_components(sir_snapshot);
if ~isequal(barrier_observed, sir_observed) || ...
        ~isequal(barrier_hidden, sir_hidden)
    error('The snapshots must use the same observed and hidden coordinates.');
end

selected_observed = barrier_observed(1:min(9, numel(barrier_observed)));
selected_hidden = barrier_hidden(1:min(6, numel(barrier_hidden)));
selected_dims = [selected_observed, selected_hidden];
if isempty(selected_dims)
    selected_dims = 1:min(15, Dx);
end

num_selected = numel(selected_dims);
num_observed = numel(selected_observed);
num_hidden = numel(selected_hidden);
num_bins = 50;

[barrier_obs_x, barrier_obs_y] = current_observation( ...
    barrier_snapshot, selected_dims, barrier_observed);
[sir_obs_x, sir_obs_y] = current_observation( ...
    sir_snapshot, selected_dims, sir_observed);

% Use one vertical range and one set of bins for both filters.
range_values = [ ...
    barrier_snapshot.particles(selected_dims, :)'; ...
    sir_snapshot.particles(selected_dims, :)'];
range_values = range_values(:);
if isfield(barrier_snapshot, 'true_state') && ...
        ~isempty(barrier_snapshot.true_state)
    range_values = [range_values; ...
        reshape(barrier_snapshot.true_state(selected_dims), [], 1)];
end
if isfield(sir_snapshot, 'true_state') && ~isempty(sir_snapshot.true_state)
    range_values = [range_values; ...
        reshape(sir_snapshot.true_state(selected_dims), [], 1)];
end
range_values = [range_values; barrier_obs_y(:); sir_obs_y(:)];

y_min = min(range_values, [], 'all');
y_max = max(range_values, [], 'all');
if y_min == y_max
    y_min = y_min - 1;
    y_max = y_max + 1;
else
    y_padding = 0.04 * (y_max - y_min);
    y_min = y_min - y_padding;
    y_max = y_max + y_padding;
end
edges = linspace(y_min, y_max, num_bins + 1);

[barrier_counts, barrier_density] = density_matrices( ...
    barrier_snapshot, selected_dims, edges);
[sir_counts, sir_density] = density_matrices( ...
    sir_snapshot, selected_dims, edges);

shared_max_count = max([barrier_counts(:); sir_counts(:)]);
if shared_max_count == 0
    shared_max_count = 1;
end
shared_max_weight = max([barrier_density(:); sir_density(:)]);
if shared_max_weight == 0
    shared_max_weight = 1;
end

% Wide, short proportions suitable for a two-column paper figure.
scrsz = get(0, 'ScreenSize');
fig_w = 1200;
fig_h = 620;
fig_pos = [(scrsz(3) - fig_w) / 2, (scrsz(4) - fig_h) / 2, ...
    fig_w, fig_h];
f = figure('Name', 'Barrier-SIR and SIR Density Comparison', ...
    'Color', 'w', 'Position', fig_pos, ...
    'PaperUnits', 'inches', 'PaperPosition', [0 0 7.2 4.6], ...
    'PaperSize', [7.2 4.6]);

layout = tiledlayout(f, 2, 1, 'TileSpacing', 'loose', ...
    'Padding', 'compact');
colormap(f, parula(256));

ax1 = nexttile(layout);
legend1 = plot_density_panel(ax1, barrier_snapshot, selected_dims, edges, ...
    barrier_counts, barrier_density, shared_max_count, ...
    shared_max_weight, barrier_obs_x, barrier_obs_y, ...
    num_observed, num_hidden, 'Barrier-SIR');

ax2 = nexttile(layout);
legend2 = plot_density_panel(ax2, sir_snapshot, selected_dims, edges, ...
    sir_counts, sir_density, shared_max_count, shared_max_weight, ...
    sir_obs_x, sir_obs_y, num_observed, num_hidden, 'SIR');

linkaxes([ax1, ax2], 'y');
ylim([ax1, ax2], [y_min, y_max]);
xlabel(layout, 'State coordinate index');
ylabel(layout, 'State value');
drawnow;
place_legend_top_left(ax1, legend1);
place_legend_top_left(ax2, legend2);

[output_dir, ~, ~] = fileparts(barrier_file);
output_base = fullfile(output_dir, sprintf( ...
    'Barrier_SIR_vs_SIR_CustomDensity_Dx%d_t%.3f', Dx, time_val));
savefig(f, [output_base '.fig']);
print(f, [output_base '.eps'], '-depsc', '-vector', '-r300');

fprintf('Density comparison saved to:\n  - %s.fig\n  - %s.eps\n', ...
    output_base, output_base);


function snapshot = load_snapshot(file_path)
if ~exist(file_path, 'file')
    error('File not found: %s', file_path);
end
data = load(file_path);
if ~isfield(data, 'snapshot')
    error('The file does not contain a snapshot structure: %s', file_path);
end
snapshot = data.snapshot;
end


function [observed, hidden] = get_components(snapshot)
Dx = snapshot.Dx;
if isfield(snapshot, 'fixed_observed_components')
    observed = snapshot.fixed_observed_components(:)';
else
    observed = 1:min(Dx, 1);
end
if isfield(snapshot, 'non_observed_components')
    hidden = snapshot.non_observed_components(:)';
else
    hidden = setdiff(1:Dx, observed);
end
end


function [obs_x, obs_y] = current_observation( ...
        snapshot, selected_dims, observed_components)
obs_x = [];
obs_y = [];
required_fields = {'ze_sparse', 'n_obs', 'step_idx'};
if ~all(isfield(snapshot, required_fields)) || isempty(snapshot.ze_sparse)
    return;
end
if snapshot.step_idx <= 0 || rem(snapshot.step_idx, snapshot.n_obs) ~= 0
    return;
end

obs_col = snapshot.step_idx / snapshot.n_obs;
if obs_col < 1 || obs_col > size(snapshot.ze_sparse, 2)
    return;
end

[is_observed, obs_row] = ismember(selected_dims, observed_components);
obs_x = find(is_observed);
if ~isempty(obs_x)
    obs_y = snapshot.ze_sparse(obs_row(is_observed), obs_col);
end
end


function [counts, density] = density_matrices( ...
        snapshot, selected_dims, edges)
num_bins = numel(edges) - 1;
num_selected = numel(selected_dims);
counts = zeros(num_bins, num_selected);
density = zeros(num_bins, num_selected);
weights = snapshot.weights(:);

for local_d = 1:num_selected
    values = snapshot.particles(selected_dims(local_d), :)';
    bin_index = discretize(values, edges);
    valid = ~isnan(bin_index);
    counts(:, local_d) = accumarray(bin_index(valid), 1, ...
        [num_bins, 1], @sum, 0);
    density(:, local_d) = accumarray(bin_index(valid), weights(valid), ...
        [num_bins, 1], @sum, 0);
end
end


function lgd = plot_density_panel(ax, snapshot, selected_dims, edges, ...
        counts, density, max_count, max_weight, obs_x, obs_y, ...
        num_observed, num_hidden, panel_title)
num_selected = numel(selected_dims);
cmap = parula(256);
lgd = gobjects(0);
hold(ax, 'on');

for local_d = 1:num_selected
    for b = 1:size(counts, 1)
        particle_count = counts(b, local_d);
        if particle_count == 0
            continue;
        end

        width = 0.8 * particle_count / max_count;
        x_left = local_d - 0.4;
        x_right = x_left + width;
        weight_sum = density(b, local_d);
        color_index = round(1 + 255 * weight_sum / max_weight);
        color_index = max(1, min(256, color_index));

        patch(ax, [x_left, x_right, x_right, x_left], ...
            [edges(b), edges(b), edges(b + 1), edges(b + 1)], ...
            cmap(color_index, :), 'EdgeColor', 'none', ...
            'FaceAlpha', 0.9, 'HandleVisibility', 'off');
    end
end

% Dotted guides pass vertically through the true/observation markers.
for local_d = 1:num_selected
    xline(ax, local_d, ':', 'Color', [0.35 0.35 0.35], ...
        'LineWidth', 0.75, 'HandleVisibility', 'off');
end

legend_handles = gobjects(0);
legend_labels = {};
if isfield(snapshot, 'true_state') && ~isempty(snapshot.true_state)
    true_values = snapshot.true_state(selected_dims);
    h_true = plot(ax, 1:num_selected, true_values, 'kx', ...
        'MarkerSize', 7, 'LineWidth', 1.4, 'LineStyle', 'none');
    legend_handles(end + 1) = h_true;
    legend_labels{end + 1} = 'True';
end
if ~isempty(obs_x)
    h_obs = plot(ax, obs_x, obs_y, 'd', 'Color', [0 0.45 0.2], ...
        'MarkerFaceColor', [0.2 0.75 0.35], 'MarkerSize', 6, ...
        'LineWidth', 1.0, 'LineStyle', 'none');
    legend_handles(end + 1) = h_obs;
    legend_labels{end + 1} = 'Observation';
end

clim(ax, [0, max_weight]);
cb = colorbar(ax);
cb.Label.String = 'Weights';

set(ax, 'XTick', 1:num_selected, ...
    'XTickLabel', arrayfun(@num2str, selected_dims, ...
    'UniformOutput', false), 'Layer', 'top');
xlim(ax, [0.5, num_selected + 0.5]);
title(ax, panel_title);
grid(ax, 'off');
box(ax, 'on');

if ~isempty(legend_handles)
    lgd = legend(ax, legend_handles, legend_labels, ...
        'Orientation', 'horizontal', 'Box', 'off', ...
        'AutoUpdate', 'off');
end

add_coordinate_group_labels(ax, num_observed, num_hidden, num_selected);
hold(ax, 'off');
end


function place_legend_top_left(ax, lgd)
if isempty(lgd) || ~isgraphics(lgd)
    return;
end
ax.Units = 'normalized';
lgd.Units = 'normalized';
legend_position = lgd.Position;
legend_position(1) = ax.Position(1);
legend_position(2) = ax.Position(2) + ax.Position(4) + 0.006;
lgd.Position = legend_position;
end


function add_coordinate_group_labels( ...
        ax, num_observed, num_hidden, num_selected)
label_y = -0.12;
if num_observed > 0
    observed_center = (num_observed / 2) / num_selected;
    text(ax, observed_center, label_y, 'observed coordinates', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontAngle', 'italic', ...
        'Clipping', 'off');
end
if num_hidden > 0
    hidden_center = (num_observed + num_hidden / 2) / num_selected;
    text(ax, hidden_center, label_y, 'hidden coordinates', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', 'FontAngle', 'italic', ...
        'Clipping', 'off');
end
end

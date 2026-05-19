% =========================================================================
% Script to analyze and visualize particle diversity from snapshots
% =========================================================================

% 1. Specify the path to your generated .mat snapshot file
file_path = 'dene/snapshots/Dx_100/Barrier_SIR_step_01000.mat'; % Change this!

% Load the data
if ~exist(file_path, 'file')
    error('File not found: %s. Please check the path.', file_path);
end
data = load(file_path);
snapshot = data.snapshot;

particles = snapshot.particles; % Dx x N
weights = snapshot.weights;     % 1 x N
Dx = snapshot.Dx;
step_idx = snapshot.step_idx;
N = length(weights);

% Simulation parameters for dynamic naming
he = 1e-3; 
time_val = step_idx * he;

% Extract Filter Name from filename
[output_dir, fname, ~] = fileparts(file_path);
parts = strsplit(fname, '_step_');
if length(parts) >= 2
    filter_name = parts{1};
else
    filter_name = 'Filter';
end
filter_title = strrep(filter_name, '_', '-'); % Avoid subscript formatting in titles

% Component Selection
if isfield(snapshot, 'fixed_observed_components')
    observed_components = snapshot.fixed_observed_components(:)';
else
    observed_components = 1:min(Dx, 1);
end
if isfield(snapshot, 'non_observed_components')
    non_observed_components = snapshot.non_observed_components(:)';
else
    non_observed_components = setdiff(1:Dx, observed_components);
end

fprintf('--- Analyzing %s Snapshot at Step %d (Time = %.3f) ---\n', filter_title, step_idx, time_val);

% =========================================================================
% PART 1: Quantitative Metrics (ESS and Variance)
% =========================================================================
ESS = 1 / sum(weights.^2);
fprintf('Effective Sample Size (ESS): %.2f (out of %d particles)\n', ESS, N);

mu = particles * weights'; 
p_diff = particles - mu; 
state_variances = sum(weights .* (p_diff.^2), 2); 
mean_system_variance = mean(state_variances);
fprintf('Mean System-wide Variance: %.4f\n\n', mean_system_variance);


% =========================================================================
% PART 2: Visualizations (Centered on Screen)
% =========================================================================
% Get screen size to calculate center positions
scrsz = get(0, 'ScreenSize'); % [left, bottom, width, height]

% -------------------------------------------------------------------------
% i) Visualization A: Bar Plot of Particle Weights
% -------------------------------------------------------------------------
fig1_w = 800; fig1_h = 400;
pos1 = [(scrsz(3)-fig1_w)/2, (scrsz(4)-fig1_h)/2, fig1_w, fig1_h];

f1 = figure('Name', 'Particle Weights', 'Color', 'w', 'Position', pos1);
bar(1:N, weights, 'FaceColor', [0.2 0.6 0.8], 'EdgeColor', 'none');
title(sprintf('%s | Dx=%d | t=%.3f\nIndividual Particle Weights', filter_title, Dx, time_val));
xlabel('Particle Index (1 to N)');
ylabel('Particle Weight');
xlim([0 N+1]);
grid on;

filename_f1 = fullfile(output_dir, sprintf('%s_WeightBarPlot_Dx%d_t%.3f.fig', filter_name, Dx, time_val));
savefig(f1, filename_f1);


% -------------------------------------------------------------------------
% ii) Visualization B: 2D Particle Spread (2 selected states)
% -------------------------------------------------------------------------
dim1 = observed_components(1);
if isempty(non_observed_components)
    dim2 = observed_components(min(2, numel(observed_components)));
else
    dim2 = non_observed_components(1);
end

fig2_w = 600; fig2_h = 500;
pos2 = [(scrsz(3)-fig2_w)/2, (scrsz(4)-fig2_h)/2, fig2_w, fig2_h];

f2 = figure('Name', '2D Particle Spread', 'Color', 'w', 'Position', pos2);
norm_w = weights / max(weights);
sizes = 10 + norm_w * 140;

scatter(particles(dim1,:), particles(dim2,:), sizes, weights, 'filled', 'MarkerFaceAlpha', 0.7);
colormap(f2, jet);
cb = colorbar;
ylabel(cb, 'Particle Weight');
title(sprintf('%s | Dx=%d | t=%.3f\nParticle Spread: X_{%d} vs X_{%d} (ESS: %.1f / %d)', ...
    filter_title, Dx, time_val, dim1, dim2, ESS, N));
xlabel(sprintf('State X_{%d}', dim1));
ylabel(sprintf('State X_{%d}', dim2));
grid on;

filename_f2 = fullfile(output_dir, sprintf('%s_Scatter2D_Dx%d_t%.3f.fig', filter_name, Dx, time_val));
savefig(f2, filename_f2);


% -------------------------------------------------------------------------
% iii) Visualization C: Custom Density Plot (Width=Count, Color=Weight)
% -------------------------------------------------------------------------
fig3_w = 1150; fig3_h = 500;
pos3 = [(scrsz(3)-fig3_w)/2, (scrsz(4)-fig3_h)/2, fig3_w, fig3_h];

f3 = figure('Name', 'Custom Density Plot', 'Color', 'w', 'Position', pos3);

sel_obs = observed_components(1:min(9, numel(observed_components)));
sel_non_obs = non_observed_components(1:min(6, numel(non_observed_components)));
selected_dims = [sel_obs, sel_non_obs];
if isempty(selected_dims)
    selected_dims = 1:min(15, Dx);
end
num_selected = numel(selected_dims);

num_bins = 50; 
min_val = min(particles(selected_dims, :), [], 'all');
max_val = max(particles(selected_dims, :), [], 'all');
if isfield(snapshot, 'true_state') && ~isempty(snapshot.true_state)
    min_val = min(min_val, min(snapshot.true_state(selected_dims)));
    max_val = max(max_val, max(snapshot.true_state(selected_dims)));
end
edges = linspace(min_val, max_val, num_bins+1);

% Matrices to store BOTH count and total weight
count_matrix = zeros(num_bins, num_selected);
density_matrix = zeros(num_bins, num_selected);

for local_d = 1:num_selected
    d = selected_dims(local_d);
    for b = 1:num_bins
        idx = particles(d,:) >= edges(b) & particles(d,:) < edges(b+1);
        count_matrix(b, local_d) = sum(idx);
        density_matrix(b, local_d) = sum(weights(idx));
    end
end

% Setup for custom patches
max_count = max(count_matrix(:)); 
if max_count == 0; max_count = 1; end % Avoid divide by zero

c_min = 0; % Minimum probability mass is 0
c_max = max(density_matrix(:));
if c_max == 0; c_max = 1; end % Avoid divide by zero for colormap scaling

cmap = colormap(f3, parula(256));
hold on;

% Draw rectangles for each bin
for local_d = 1:num_selected
    for b = 1:num_bins
        c = count_matrix(b, local_d);
        if c > 0
            w_sum = density_matrix(b, local_d);
            
            % Width proportional to number of particles (Max width = 0.8 to prevent overlap)
            width = 0.8 * (c / max_count);
            x_left = local_d - width/2;
            x_right = local_d + width/2;
            y_bottom = edges(b);
            y_top = edges(b+1);
            
            % Map weight sum to colormap index (1 to 256)
            c_idx = round(1 + 255 * (w_sum - c_min) / (c_max - c_min));
            c_idx = max(1, min(256, c_idx));
            face_color = cmap(c_idx, :);
            
            patch([x_left, x_right, x_right, x_left], ...
                  [y_bottom, y_bottom, y_top, y_top], ...
                  face_color, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
        end
    end
end

% Set up colorbar for the patches
clim([c_min c_max]);
cb = colorbar;
ylabel(cb, 'Probability Mass (Sum of Weights)');

% Overlays (Weighted mean, True State, Observed Value) - REMOVED DisplayNames for Legend
mu_selected = mu(selected_dims);
plot(1:num_selected, mu_selected, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);

if isfield(snapshot, 'true_state') && ~isempty(snapshot.true_state)
    true_selected = snapshot.true_state(selected_dims);
    plot(1:num_selected, true_selected, 'kx', 'MarkerSize', 7, 'LineWidth', 1.5);
end

if isfield(snapshot, 'ze_sparse') && isfield(snapshot, 'n_obs') && ~isempty(snapshot.ze_sparse)
    obs_col = floor(step_idx / snapshot.n_obs);
    if obs_col >= 1 && obs_col <= size(snapshot.ze_sparse, 2)
        [is_obs, obs_local_idx] = ismember(selected_dims, observed_components);
        obs_x = find(is_obs);
        if ~isempty(obs_x)
            obs_y = snapshot.ze_sparse(obs_local_idx(is_obs), obs_col);
            plot(obs_x, obs_y, 'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
        end
    end
end
hold off;

% Formatting
y_tick_labels = arrayfun(@(d) sprintf('%d', d), selected_dims, 'UniformOutput', false);
set(gca, 'XTick', 1:num_selected, 'XTickLabel', y_tick_labels);
title(sprintf('%s | Dx=%d | t=%.3f\nCustom Density: Width = Particle Count, Color = Probability Mass', ...
    filter_title, Dx, time_val));
xlabel('Selected State Dimension Index');
ylabel('State Value');
xlim([0.5 num_selected+0.5]);
ylim([min_val max_val]);
grid on;

% Save
filename_f3 = fullfile(output_dir, sprintf('%s_CustomDensity_Dx%d_t%.3f.fig', filter_name, Dx, time_val));
savefig(f3, filename_f3);

fprintf('Plots successfully saved as .fig files to:\n  - %s\n  - %s\n  - %s\n', filename_f1, filename_f2, filename_f3);
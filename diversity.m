% =========================================================================
% Script to analyze and visualize particle diversity from snapshots
% =========================================================================

% 1. Specify the path to your generated .mat snapshot file
%file_path = '/snapshots/Dx_500/Barrier_SIR_step_01400.mat'; % Change this!
file_path = 'small/snapshots/Dx_10/Barrier_SIR_step_00010.mat';
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

fprintf('--- Analyzing Snapshot at Step %d ---\n', step_idx);

% =========================================================================
% PART 1: Quantitative Metrics (ESS and Variance)
% =========================================================================

% A. Effective Sample Size (ESS)
ESS = 1 / sum(weights.^2);
fprintf('Effective Sample Size (ESS): %.2f (out of %d particles)\n', ESS, N);

% B. Spatial Spread (Weighted Variance)
% Weighted mean (Dx x 1)
mu = particles * weights'; 

% Difference from mean (Dx x N)
p_diff = particles - mu; 

% Variance for each state dimension (Dx x 1)
state_variances = sum(weights .* (p_diff.^2), 2); 

% Average variance across all dimensions
mean_system_variance = mean(state_variances);
fprintf('Mean System-wide Variance: %.4f\n\n', mean_system_variance);


% =========================================================================
% PART 2: Visualizations
% =========================================================================

% --- Visualization A: 2D State Scatter Plot ---
% Use first observed component and first non-observed component
dim1 = observed_components(1);
if isempty(non_observed_components)
    dim2 = observed_components(min(2, numel(observed_components)));
else
    dim2 = non_observed_components(1);
end

figure('Name', '2D Particle Spread', 'Color', 'w', 'Position', [100, 100, 600, 500]);
% Normalize weights to scale marker sizes (min size 10, max size 150)
norm_w = weights / max(weights);
sizes = 10 + norm_w * 140;

% Scatter plot colored and sized by weight
scatter(particles(dim1,:), particles(dim2,:), sizes, weights, 'filled', 'MarkerFaceAlpha', 0.7);
colormap(jet);
cb = colorbar;
ylabel(cb, 'Particle Weight');
title(sprintf('Particle Spread: X_{%d} vs X_{%d} (Step %d)\nESS: %.1f / %d', dim1, dim2, step_idx, ESS, N));
xlabel(sprintf('State X_{%d}', dim1));
ylabel(sprintf('State X_{%d}', dim2));
grid on;


% --- Visualization B: 2D Marginal Density (Weighted Histogram across selected states) ---
figure('Name', 'Weighted Histogram All States', 'Color', 'w', 'Position', [750, 100, 1000, 500]);

% Keep up to 15 components: 9 observed + 6 non-observed
sel_obs = observed_components(1:min(9, numel(observed_components)));
sel_non_obs = non_observed_components(1:min(6, numel(non_observed_components)));
selected_dims = [sel_obs, sel_non_obs];
if isempty(selected_dims)
    selected_dims = 1:min(15, Dx);
end
num_selected = numel(selected_dims);

% Define common value bins for all states based on the global min and max
num_bins = 50; 
min_val = min(particles(selected_dims, :), [], 'all');
max_val = max(particles(selected_dims, :), [], 'all');
if isfield(snapshot, 'true_state') && ~isempty(snapshot.true_state)
    min_val = min(min_val, min(snapshot.true_state(selected_dims)));
    max_val = max(max_val, max(snapshot.true_state(selected_dims)));
end
edges = linspace(min_val, max_val, num_bins);
centers = (edges(1:end-1) + edges(2:end)) / 2;

% Preallocate density matrix (rows = value bins, cols = state dimensions)
density_matrix = zeros(length(centers), num_selected);

% Compute weighted histogram for each state dimension
for local_d = 1:num_selected
    d = selected_dims(local_d);
    for b = 1:length(centers)
        % Find particles whose d-th state falls in the b-th bin
        idx = particles(d,:) >= edges(b) & particles(d,:) < edges(b+1);
        density_matrix(b, local_d) = sum(weights(idx));
    end
end

% Plot as a heatmap
imagesc(1:num_selected, centers, density_matrix);
set(gca, 'YDir', 'normal'); % Keep Y-axis increasing upwards
colormap(parula); % 'jet' or 'hot' also look good for density
cb = colorbar;
ylabel(cb, 'Probability Mass (Sum of Weights)');

% Overlay the weighted mean for reference (Points only, no connecting lines)
hold on;
mu_selected = mu(selected_dims);
plot(1:num_selected, mu_selected, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Weighted Mean');

% Overlay true state values for exact signal
if isfield(snapshot, 'true_state') && ~isempty(snapshot.true_state)
    true_selected = snapshot.true_state(selected_dims);
    plot(1:num_selected, true_selected, 'kx', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'True State');
end

% For observed components, also overlay the observed value at this step
if isfield(snapshot, 'ze_sparse') && isfield(snapshot, 'n_obs') && ~isempty(snapshot.ze_sparse)
    obs_col = floor(step_idx / snapshot.n_obs);
    if obs_col >= 1 && obs_col <= size(snapshot.ze_sparse, 2)
        [is_obs, obs_local_idx] = ismember(selected_dims, observed_components);
        obs_x = find(is_obs);
        if ~isempty(obs_x)
            obs_y = snapshot.ze_sparse(obs_local_idx(is_obs), obs_col);
            plot(obs_x, obs_y, 'gd', 'MarkerFaceColor', 'g', 'MarkerSize', 6, 'DisplayName', 'Observed Value');
        end
    end
end
hold off;
y_tick_labels = arrayfun(@(d) sprintf('%d', d), selected_dims, 'UniformOutput', false);
set(gca, 'XTick', 1:num_selected, 'XTickLabel', y_tick_labels);
title(sprintf('Weighted Histogram Density Across %d Selected States (Step %d)', num_selected, step_idx));
xlabel('Selected State Dimension Index');
ylabel('State Value');
xlim([1 num_selected]);
ylim([min_val max_val]);
legend('Location', 'northeast');

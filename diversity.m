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
% Choose two adjacent dimensions (e.g., X1 and X2)
dim1 = 1; 
dim2 = 2;

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


% --- Visualization B: 2D Marginal Density (Weighted Histogram across all states) ---
figure('Name', 'Weighted Histogram All States', 'Color', 'w', 'Position', [750, 100, 1000, 500]);

% Define common value bins for all states based on the global min and max
num_bins = 50; 
min_val = min(particles(:));
max_val = max(particles(:));
edges = linspace(min_val, max_val, num_bins);
centers = (edges(1:end-1) + edges(2:end)) / 2;

% Preallocate density matrix (rows = value bins, cols = state dimensions)
density_matrix = zeros(length(centers), Dx);

% Compute weighted histogram for each state dimension
for d = 1:Dx
    for b = 1:length(centers)
        % Find particles whose d-th state falls in the b-th bin
        idx = particles(d,:) >= edges(b) & particles(d,:) < edges(b+1);
        density_matrix(b, d) = sum(weights(idx));
    end
end

% Plot as a heatmap
imagesc(1:Dx, centers, density_matrix);
set(gca, 'YDir', 'normal'); % Keep Y-axis increasing upwards
colormap(parula); % 'jet' or 'hot' also look good for density
cb = colorbar;
ylabel(cb, 'Probability Mass (Sum of Weights)');

% Overlay the weighted mean for reference (Points only, no connecting lines)
hold on;
plot(1:Dx, mu, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Weighted Mean');
hold off;
title(sprintf('Weighted Histogram Density Across All %d States (Step %d)', Dx, step_idx));
xlabel('State Dimension Index (1 to Dx)');
ylabel('State Value');
xlim([1 Dx]);
ylim([min_val max_val]);
legend('Location', 'northeast');
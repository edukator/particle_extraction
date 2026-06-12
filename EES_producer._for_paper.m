% EES_PRODUCER Average final-time effective sample size across sim1,...,sim10.
%
% Reads:
%   simN/snapshots/Dx_D/Barrier_SIR_step_10000.mat
%   simN/snapshots/Dx_D/SIR_step_10000.mat
%
% Produces:
%   EES_results.txt
%   Barrier_SIR_average_ESS_vs_Dx.fig
%   Barrier_SIR_average_ESS_vs_Dx.png
%   SIR_average_ESS_vs_Dx.fig
%   SIR_average_ESS_vs_Dx.png

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

simulation_ids = 1:10;
dx_values = [10, 50, 100, 250, 500, 750, 1000, 1250];
filter_file_names = {'Barrier_SIR', 'SIR'};
filter_display_names = {'Barrier-SIR', 'SIR'};
final_step = 10000;

num_simulations = numel(simulation_ids);
num_dx = numel(dx_values);
num_filters = numel(filter_file_names);
ess_values = nan(num_simulations, num_dx, num_filters);

for sim_idx = 1:num_simulations
    sim_id = simulation_ids(sim_idx);

    for dx_idx = 1:num_dx
        Dx = dx_values(dx_idx);

        for filter_idx = 1:num_filters
            snapshot_path = fullfile( ...
                script_dir, sprintf('sim%d', sim_id), 'snapshots', ...
                sprintf('Dx_%d', Dx), ...
                sprintf('%s_step_%05d.mat', ...
                    filter_file_names{filter_idx}, final_step));

            if ~isfile(snapshot_path)
                error('EESProducer:MissingFile', ...
                    'Snapshot file not found: %s', snapshot_path);
            end

            loaded_data = load(snapshot_path, 'snapshot');
            if ~isfield(loaded_data, 'snapshot') || ...
                    ~isfield(loaded_data.snapshot, 'weights')
                error('EESProducer:MissingWeights', ...
                    'Variable snapshot.weights is missing from: %s', ...
                    snapshot_path);
            end

            snapshot = loaded_data.snapshot;
            weights = double(snapshot.weights(:));

            if isempty(weights) || any(~isfinite(weights)) || any(weights < 0)
                error('EESProducer:InvalidWeights', ...
                    'Invalid particle weights in: %s', snapshot_path);
            end

            weight_sum = sum(weights);
            if weight_sum <= 0
                error('EESProducer:InvalidWeights', ...
                    'Particle weights sum to zero in: %s', snapshot_path);
            end

            % Saved filter weights should already sum to one. Normalizing here
            % avoids small floating-point deviations affecting the ESS.
            weights = weights ./ weight_sum;
            ess_values(sim_idx, dx_idx, filter_idx) = ...
                1 / sum(weights .^ 2);
        end
    end
end

average_ess = reshape(mean(ess_values, 1), num_dx, num_filters);
std_ess = reshape(std(ess_values, 0, 1), num_dx, num_filters);
barrier_average_ess = average_ess(:, 1);
barrier_std_ess = std_ess(:, 1);
sir_average_ess = average_ess(:, 2);
sir_std_ess = std_ess(:, 2);

barrier_table = table(dx_values(:), barrier_average_ess, barrier_std_ess, ...
    'VariableNames', {'Dx', 'Average_ESS', 'Std_ESS'});
sir_table = table(dx_values(:), sir_average_ess, sir_std_ess, ...
    'VariableNames', {'Dx', 'Average_ESS', 'Std_ESS'});
combined_table = table( ...
    dx_values(:), barrier_average_ess, barrier_std_ess, ...
    sir_average_ess, sir_std_ess, ...
    'VariableNames', {'Dx', 'Barrier_SIR_Average_ESS', ...
    'Barrier_SIR_Std_ESS', 'SIR_Average_ESS', 'SIR_Std_ESS'});

fprintf('\nBarrier-SIR: Dx versus average final-time ESS (%d simulations)\n', ...
    num_simulations);
disp(barrier_table);

fprintf('SIR: Dx versus average final-time ESS (%d simulations)\n', ...
    num_simulations);
disp(sir_table);

fprintf('Combined average final-time ESS table\n');
disp(combined_table);

output_path = fullfile(script_dir, 'EES_results.txt');
output_file = fopen(output_path, 'w');
if output_file == -1
    error('EESProducer:OutputFile', ...
        'Could not open output file for writing: %s', output_path);
end
output_cleanup = onCleanup(@() fclose(output_file));

fprintf(output_file, ...
    'Average final-time effective sample size (ESS) over sim1 to sim10\n');
fprintf(output_file, 'Snapshot step: %d\n\n', final_step);

write_ess_table(output_file, ...
    'Barrier-SIR: Dx versus average ESS', dx_values, ...
    barrier_average_ess, barrier_std_ess);
write_ess_table(output_file, ...
    'SIR: Dx versus average ESS', dx_values, ...
    sir_average_ess, sir_std_ess);

fprintf(output_file, 'Combined table\n');
fprintf(output_file, '%-8s %-26s %-22s %-20s %-16s\n', ...
    'Dx', 'Barrier_SIR_Average_ESS', 'Barrier_SIR_Std_ESS', ...
    'SIR_Average_ESS', 'SIR_Std_ESS');
fprintf(output_file, '%s\n', repmat('-', 1, 96));
for dx_idx = 1:num_dx
    fprintf(output_file, '%-8d %-26.6f %-22.6f %-20.6f %-16.6f\n', ...
        dx_values(dx_idx), barrier_average_ess(dx_idx), ...
        barrier_std_ess(dx_idx), sir_average_ess(dx_idx), ...
        sir_std_ess(dx_idx));
end

clear output_cleanup;

create_ess_plot(dx_values, barrier_average_ess, ...
    filter_display_names{1}, ...
    fullfile(script_dir, 'Barrier_SIR_average_ESS_vs_Dx'));
create_ess_plot(dx_values, sir_average_ess, ...
    filter_display_names{2}, ...
    fullfile(script_dir, 'SIR_average_ESS_vs_Dx'));

fprintf('ESS results written to: %s\n', output_path);
fprintf('Barrier-SIR and SIR figures saved as .fig and .png files.\n');

function write_ess_table( ...
        file_id, title_text, dx_values, average_values, std_values)
    fprintf(file_id, '%s\n', title_text);
    fprintf(file_id, '%-8s %-20s %-20s\n', ...
        'Dx', 'Average_ESS', 'Std_ESS');
    fprintf(file_id, '%s\n', repmat('-', 1, 50));
    for idx = 1:numel(dx_values)
        fprintf(file_id, '%-8d %-20.6f %-20.6f\n', ...
            dx_values(idx), average_values(idx), std_values(idx));
    end
    fprintf(file_id, '\n');
end

function create_ess_plot(dx_values, average_values, filter_name, output_stem)
    figure_handle = figure('Color', 'w', 'Visible', 'off');
    dx_positions = 1:numel(dx_values);
    plot(dx_positions, average_values, '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 7, ...
        'MarkerFaceColor', [0.2, 0.45, 0.75]);
    xlabel('Dx');
    ylabel('Average ESS');
    title(sprintf('%s: Average Final-Time ESS vs Dx', filter_name));
    grid on;
    xticks(dx_positions);
    xticklabels(string(dx_values));
    xlim([0.75, numel(dx_values) + 0.25]);

    savefig(figure_handle, [output_stem, '.fig']);
    exportgraphics(figure_handle, [output_stem, '.png'], 'Resolution', 300);
    close(figure_handle);
end

function [average_ess_table, snapshot_count_table] = compute_average_ess(experiment_dir)
%COMPUTE_AVERAGE_ESS Average ESS over all stored snapshot times.
%
%   average_ess_table = compute_average_ess(experiment_dir)
%   [average_ess_table, snapshot_count_table] = ...
%       compute_average_ess(experiment_dir)
%
% experiment_dir can be either a main_cli output directory, such as
% "taming_experiment", or its "snapshots" subdirectory. The output table
% has one row per filter and one column per state dimension.

    if nargin < 1 || isempty(experiment_dir)
        experiment_dir = 'taming_experiment';
    end

    experiment_dir = char(experiment_dir);
    nested_snapshot_dir = fullfile(experiment_dir, 'snapshots');
    if isfolder(nested_snapshot_dir)
        snapshot_dir = nested_snapshot_dir;
    elseif isfolder(experiment_dir)
        snapshot_dir = experiment_dir;
    else
        error('compute_average_ess:DirectoryNotFound', ...
            'Snapshot directory not found: %s', experiment_dir);
    end

    dx_directories = dir(fullfile(snapshot_dir, 'Dx_*'));
    dx_directories = dx_directories([dx_directories.isdir]);

    records = struct('path', {}, 'filter_name', {}, 'Dx', {});
    for dx_idx = 1:numel(dx_directories)
        dx_token = regexp(dx_directories(dx_idx).name, ...
            '^Dx_(\d+)$', 'tokens', 'once');
        if isempty(dx_token)
            continue;
        end

        Dx = str2double(dx_token{1});
        dx_path = fullfile(dx_directories(dx_idx).folder, ...
            dx_directories(dx_idx).name);
        snapshot_files = dir(fullfile(dx_path, '*_step_*.mat'));

        for file_idx = 1:numel(snapshot_files)
            file_token = regexp(snapshot_files(file_idx).name, ...
                '^(.*)_step_\d+\.mat$', 'tokens', 'once');
            if isempty(file_token)
                continue;
            end

            record_idx = numel(records) + 1;
            records(record_idx).path = fullfile( ...
                snapshot_files(file_idx).folder, ...
                snapshot_files(file_idx).name);
            records(record_idx).filter_name = file_token{1};
            records(record_idx).Dx = Dx;
        end
    end

    if isempty(records)
        error('compute_average_ess:NoSnapshots', ...
            'No files matching *_step_*.mat were found under: %s', ...
            snapshot_dir);
    end

    dx_values = sort(unique([records.Dx]));
    filter_names = sort(unique(string({records.filter_name})));
    ess_sum = zeros(numel(filter_names), numel(dx_values));
    snapshot_count = zeros(size(ess_sum));

    for record_idx = 1:numel(records)
        loaded_data = load(records(record_idx).path, 'snapshot');
        if ~isfield(loaded_data, 'snapshot') || ...
                ~isfield(loaded_data.snapshot, 'weights')
            error('compute_average_ess:MissingWeights', ...
                'snapshot.weights is missing from: %s', ...
                records(record_idx).path);
        end

        weights = double(loaded_data.snapshot.weights(:));
        if isempty(weights) || any(~isfinite(weights)) || any(weights < 0)
            error('compute_average_ess:InvalidWeights', ...
                'Invalid particle weights in: %s', records(record_idx).path);
        end

        weight_sum = sum(weights);
        if weight_sum <= 0
            error('compute_average_ess:InvalidWeights', ...
                'Particle weights sum to zero in: %s', ...
                records(record_idx).path);
        end

        weights = weights ./ weight_sum;
        ess = 1 / sum(weights .^ 2);

        filter_idx = find(filter_names == records(record_idx).filter_name, 1);
        dx_idx = find(dx_values == records(record_idx).Dx, 1);
        ess_sum(filter_idx, dx_idx) = ess_sum(filter_idx, dx_idx) + ess;
        snapshot_count(filter_idx, dx_idx) = ...
            snapshot_count(filter_idx, dx_idx) + 1;
    end

    average_ess = ess_sum ./ snapshot_count;
    average_ess(snapshot_count == 0) = NaN;

    variable_names = cellstr(compose('Dx_%d', dx_values));
    row_names = cellstr(filter_names);
    average_ess_table = array2table(average_ess, ...
        'VariableNames', variable_names, 'RowNames', row_names);
    snapshot_count_table = array2table(snapshot_count, ...
        'VariableNames', variable_names, 'RowNames', row_names);
    average_ess_table.Properties.DimensionNames{1} = 'Filter';
    snapshot_count_table.Properties.DimensionNames{1} = 'Filter';

    fprintf('\nAverage ESS across stored snapshot times\n');
    fprintf('Snapshot directory: %s\n', snapshot_dir);
    disp(average_ess_table);
end

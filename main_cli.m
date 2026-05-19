function main_cli(sx, sz, is_N_fixed, fig_dir)
% main(sx, sz, is_N_fixed, fig_dir)
%   sx, sz      : signal and observation noise std
%   is_N_fixed  : logical, if true N is fixed, else N = DxArray
%   fig_dir     : directory where the .fig file will be saved

    % Optional but often useful: ensure output folder exists
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    %% Parameters

    % LORENZ 96 parameters
    F       = 8;        % forcing parameter
    he      = 1e-3;     % time step, standard Euler
    t_final = 10;       % duration of the simulation in natural time units
    NTe     = fix(t_final/he);  % no. of discrete time steps, standard Euler

    % sz, sx are now taken from input arguments
    % sz = ...;
    % sx = ...;

    tobs  = 0.1;                % continuous time between observations

    % Particle filter parameters
    ness_thr = 0.7;             % NESS threshold for resampling
    eps0     = 0.05;            % parameter of the hut distribution (unused below?)

    % Simulation parameters
    n_steps = ceil(5/he);       % number of time step to skip transient solution      
    n_obs   = ceil(tobs/he);    % number of subintervals between two observations

    filtered_solution_indices = 1:n_obs:NTe+1; % filtered solution is computed at these indices
    coarse_time_mesh          = he.*(0:n_obs:NTe); % time mesh corresponding to those indices

    full_indices   = 1:NTe+1;   % ground truth signal and predicted solution are computed at these indices
    fine_time_mesh = he*(0:NTe);

    % Parameters for hypercube and tube preservation passed to sir_rejectproject
    parhyper.alpha = 0.5;
    parhyper.beta  = 1;
    parhyper.Lm    = 0.8;

    %% barrier parameters
    r_obs                  = 4;

    barrier_params.mu      = 80;     
    barrier_params.p       = r_obs;
    barrier_params.k       = 100;

    dummy_params = struct();

    
     DxArray=[100];
      save_step_indices = [10,1000];
    snapshot_root_dir = fullfile(fig_dir, 'snapshots');
    if ~exist(snapshot_root_dir, 'dir')
        mkdir(snapshot_root_dir);
    end
 

    filters = {
        SIRFilter(dummy_params,     DxArray), ...
        ENKFFilter(dummy_params,    DxArray), ...
        ENKF_Barrier_Filter(dummy_params, DxArray),	...
        APFFilter(dummy_params,     DxArray), ...
        BarrierFilter(dummy_params, DxArray)
    };

    % NArray: fixed or equal to DxArray, depending on input flag
    if is_N_fixed
        NArray = 500*ones(size(DxArray));
    else
        NArray = DxArray;
    end

    %% Main loop over dimensions
    for iDX = 1:numel(DxArray)
        Dx = DxArray(iDX);
        N  = NArray(iDX);

        parhyper.r0 = parhyper.Lm^(1/Dx);

        Dz = fix(3*Dx/5);          % set full observation
        % Dz = Dx;

       

        fixed_observed_components = randsample(Dx, Dz);
        fixed_observed_components = sort(fixed_observed_components);
            ok = 0;
            while ~ok
                n_steps = ceil(5/he);           % up to half of final time ?
                Wx0     = sqrt(he)*randn([Dx n_steps]); % Brownian increment

                [x_ini, ~] = exp_euler(rand([Dx 1]), he, F, n_steps, Dx, Wx0, sx);

                idx = randsample(fix(n_steps/2):n_steps, 1); % choose an integer in [n_steps/2, n_steps]
                x0  = x_ini(:, idx);

                % now we run with the 'regular' initialisation
                Wx = sqrt(he)*randn([Dx NTe]);  % create new increment
                [x, ok] = exp_euler(x0, he, F, NTe, Dx, Wx, sx);
                % if ok == 1, x does not contain NaN values
            end

            H0   = eye(Dx) + randn([Dx Dx]).*5e-4;
            H0x  = H0 * x(1:Dx, (n_obs+1):n_obs:NTe+1); % skip initial time, start from t=h
            ze_full = H0x + sz*randn(size(H0x));

            H         = H0(fixed_observed_components, :);       % observation matrix
            ze_sparse = ze_full(fixed_observed_components, :);  % observed values

            X0 = x0 + sx*randn([Dx N]);  % initial particles

            truth      = x(:, filtered_solution_indices);
            truth_full = x;

            Pd_f = mean(sum(x(:, filtered_solution_indices).^2));
            Pd_p = mean(sum(x.^2));

            params = struct( ...
                'F',                        F, ...
                'Dx',                       Dx, ...
                'sx',                       sx, ...
                'sz',                       sz, ...
                'he',                       he, ...
                'NTe',                      NTe, ...
                'n_obs',                    n_obs, ...
                'ze_sparse',                ze_sparse, ...
                'H',                        H, ...
                'X0',                       X0, ...
                'ness_thr',                 ness_thr, ...
                'Dz',                       Dz, ...
                'fixed_observed_components', fixed_observed_components, ...
                'truth',                    truth, ...
                'truth_full',               truth_full, ...
                'Pd_f',                     Pd_f, ...
                'Pd_p',                     Pd_p, ...
                'parhyper',                 parhyper, ...
                'barrier_params',           barrier_params, ...
                'save_step_indices',         save_step_indices, ...
                'snapshot_dir',             fullfile(snapshot_root_dir, sprintf('Dx_%d', Dx)) ...
            );

            % Assign updated params to filters and run them
            for k = 1:numel(filters)
                filters{k}.params = params;
                filters{k}.run(iDX);
            end
    end % Dx loop

    % Display metrics
    for k = 1:numel(filters)
        disp(filters{k}.metric);
    end

    fprintf('Particle/weight snapshots are saved on the fly to: %s\n', snapshot_root_dir);

    %% Plot and save as .fig (invisible)
    f = figure('Visible', 'on');
    hold on;

    for k = 1:numel(filters)
        metric = filters{k}.metric;
        dx_values = metric.DxArray(:)';
        mean_mse = metric.MSEf(:)';
        std_mse  = zeros(size(mean_mse));

        line_handle = plot(dx_values, mean_mse, '-o', 'DisplayName', metric.name);

        upper = mean_mse + std_mse;
        lower = mean_mse - std_mse;

        band_handle = fill([dx_values, fliplr(dx_values)], [upper, fliplr(lower)], ...
            get(line_handle, 'Color'), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        set(band_handle, 'HandleVisibility', 'off');
        uistack(band_handle, 'bottom');
    end

    xlabel('Dx');
    ylabel('Mean MSE_f');
    title('Mean Filtered MSE vs Dx');
    legend('Location', 'best');
    grid on;

    fig_path = fullfile(fig_dir, 'mean_filtered_MSE_vs_Dx.fig');
    savefig(f, fig_path);
    close(f);

    fprintf('Figure saved to: %s\n', fig_path);
end

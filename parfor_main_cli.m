function parfor_main_cli(sx, sz, NR, is_N_fixed, fig_dir, in_cluster)
% parfor_main_cli(sx, sz, NR, is_N_fixed, fig_dir, in_cluster)
%   Parallel version using parfor over simulation runs (nr) for each Dx.
%   sx, sz      : signal and observation noise std
%   NR          : number of runs per Dx
%   is_N_fixed  : logical, if true N is fixed, else N = DxArray
%   fig_dir     : directory where the .fig file will be saved
%   in_cluster  : logical, if true do not show plot; only save

    % Optional but often useful: ensure output folder exists
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end

    %% Parameters

    % LORENZ 96 parameters
    F       = 8;        % forcing parameter
    he      = 1e-3;     % time step, standard Euler
    t_final = 10;        % duration of the simulationufnad simulation in natural time units
    NTe     = fix(t_final/he);  % no. of discrete time steps, standard Euler

    tobs  = 0.1;                % continuous time between observations

    % Particle filter parameters
    ness_thr = 0.7;             % NESS threshold for resampling
    eps0     = 0.05;            % parameter of the hut distribution (unused below?)

    % Simulation parameters
    n_obs   = ceil(tobs/he);    % number of subintervals between two observations

    filtered_solution_indices = 1:n_obs:NTe+1; % filtered solution is computed at these indices

    % Parameters for hypercube and tube preservation passed to sir_rejectproject
    parhyper.alpha = 0.5;
    parhyper.beta  = 1;
    parhyper.Lm    = 0.8;

    %% barrier parameters
    r_obs                  = 4; % not used ?

    barrier_params.mu      = 50;
    barrier_params.p       = r_obs;
    barrier_params.k       = 100;

    dummy_params = struct();

    DxArray =  [10,50,100,250,500,750,1000,1250];
   % DxArray =  [10,50,100];
    % NArray: fixed or equal to DxArray, depending on input flag
    if is_N_fixed
        NArray = 400*ones(size(DxArray));
    else
        NArray = DxArray;
    end

    % Filter definitions
    filter_classes = {@SIRFilter, @ENKFFilter, @ENKF_Barrier_Filter, @APFFilter, @BarrierFilter,...
        @GUIDED_GIRSANOV_a0,@GUIDED_GIRSANOV_BARRIER_a0};
    filter_names   = {'SIR', 'EnKF', 'Barrier-ENKF', 'APF', 'Barrier-SIR',...
        'Guided','Barrier-Guided'};
    nFilters = numel(filter_classes);

    % Metrics storage
    metrics = cell(1, nFilters);
    for k = 1:nFilters
        metrics{k} = ErrorMetric(filter_names{k}, DxArray, NR);
    end

    % Preallocate result arrays for parfor assignments
    mse_f_all = nan(nFilters, numel(DxArray), NR);
    mse_p_all = nan(nFilters, numel(DxArray), NR);
    runtime_all = nan(nFilters, numel(DxArray), NR);

    %% Main loop over dimensions and runs
    for iDX = 1:numel(DxArray)
        Dx = DxArray(iDX);
        N  = NArray(iDX);

        parhyper_current = parhyper;
        parhyper_current.r0 = parhyper.Lm^(1/Dx);

        Dz = fix(3*Dx/5);          % set full observation

        parfor nr = 1:NR
            % Pre-initialize variables flagged as temporaries to satisfy parfor analysis
            x0 = [];
            x  = [];
            n_steps = ceil(5/he);           % up to half of final time ?

            fixed_observed_components = randsample(Dx, Dz);
            fixed_observed_components = sort(fixed_observed_components);

            ok = 0;
            while ~ok
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
                'parhyper',                 parhyper_current, ...
                'barrier_params',           barrier_params ...
            );

            for k = 1:nFilters
                filter_obj = filter_classes{k}(dummy_params, DxArray, NR);
                filter_obj.params = params;
                filter_obj.run(iDX, nr);
                mse_f_all(k, iDX, nr) = filter_obj.metric.MSEf(iDX, nr);
                mse_p_all(k, iDX, nr) = filter_obj.metric.MSEp(iDX, nr);
                runtime_all(k, iDX, nr) = filter_obj.metric.runtime(iDX, nr);
            end
        end % nr parfor
    end % Dx loop

    % Consolidate results into metrics
    for k = 1:nFilters
        metrics{k}.MSEf = squeeze(mse_f_all(k, :, :));
        metrics{k}.MSEp = squeeze(mse_p_all(k, :, :));
        metrics{k}.runtime = squeeze(runtime_all(k, :, :));
    end

    % Display metrics
    for k = 1:nFilters
        disp(metrics{k});
    end

    %% Plot and save
    if in_cluster
        figVisible = 'off';
    else
        figVisible = 'on';
    end

    f = figure('Visible', figVisible);
    hold on;

   % ---- STYLE LIBRARY (>= 12 distinct styles) ----
% MATLAB default color order is good and readable; we extend it with extra colors.
colors = [ ...
    0.0000 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange
    0.9290 0.6940 0.1250;  % yellow
    0.4940 0.1840 0.5560;  % purple
    0.4660 0.6740 0.1880;  % green
    0.3010 0.7450 0.9330;  % cyan
    0.6350 0.0780 0.1840;  % dark red
    0.2500 0.2500 0.2500;  % dark gray
    0.7500 0.0000 0.7500;  % magenta-ish
    0.0000 0.6000 0.0000;  % deep green
    0.0000 0.0000 0.0000;  % black
    0.9000 0.4000 0.7000;  % pink-ish
];

markers = {'o','s','^','d','v','>','<','p','h','x','+','*'};

nStyles = min(size(colors,1), numel(markers));  % at least 12 here

% (Optional) adjust global marker/line sizing
lineWidth  = 1.5;
markerSize = 7;

for k = 1:nFilters
    metric = metrics{k};
    dx_values = metric.DxArray(:)';

    mean_mse = mean(metric.MSEf, 2, 'omitnan')';
    std_mse  = std(metric.MSEf, 0, 2, 'omitnan')';

    % pick style "in order"; if you ever have > nStyles filters, it cycles
    sidx = mod(k-1, nStyles) + 1;
    c = colors(sidx,:);
    m = markers{sidx};

    % line + markers
    line_handle = plot(dx_values, mean_mse, ...
        '-', ...
        'Color', c, ...
        'Marker', m, ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'MarkerFaceColor', c, ...
        'DisplayName', metric.name);

    % uncertainty band (same color, transparent)
    upper = mean_mse + std_mse;
    lower = mean_mse - std_mse;

    band_handle = fill([dx_values, fliplr(dx_values)], [upper, fliplr(lower)], ...
        c, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    set(band_handle, 'HandleVisibility', 'off');
    uistack(band_handle, 'bottom');
end

xlabel('d_x');
ylabel('NMSE');
%title('Mean Filtered MSE vs Dx');
legend('Location', 'best');
grid on

    fig_path = fullfile(fig_dir, 'mean_filtered_MSE_vs_Dx.fig');
    savefig(f, fig_path);

      if in_cluster
        close(f); % do not show in cluster mode
      else
         drawnow;  % ensure it renders when running interactively
      end

    fprintf('Figure saved to: %s\n', fig_path);
end

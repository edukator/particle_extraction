

clear;
clc;
close all;




	  
%% Parameters

% LORENZ 96 parameters
F = 8;			% forcing parameter
he = 1e-3;  		% time step, standard Euler
t_final = 10;        %30   % duration of the simulation in natural time units
NTe = fix(t_final/he);  % no. of discrete time steps, standard Euler

sz = 1; %sqrt(1/4)        % std of the observations; full filtering
sx = 1;         %5  ,sqrt(1/2)      % std of the signal noise / diffusion coefficient
tobs = 0.1;            % continuous time between observations
     	           
% Particle filter parameters

ness_thr = 0.7;		% NESS threshold for resampling
%
eps0 = 0.05;            % parameter of the hut distribution
%
% Simulation parameters
%
n_steps = ceil(5/he);    % number of time step to skip transient solution      
n_obs = ceil(tobs/he);   % number of subintervals between two observations
filtered_solution_indices=1:n_obs:NTe+1; % filtered solution is computed at these index 
coarse_time_mesh=he.*(0:n_obs:NTe); % time mesh corresponding those inex
is_N_fixed=true;
full_indices=1:NTe+1;    % ground truth signal and predicted solution are computed at these index
fine_time_mesh=he*(0:NTe);

% Parameters for hypercube and tube preservation passed to sir_rejectproject
% parhyper : hypercube parameters
% parmeth : method and sampling parameters 

parhyper.alpha=0.5;
parhyper.beta=1;
parhyper.Lm=0.8;




%		    %
% 
%% barrier paramemters; 
r_obs=3;
barrier_params.alpha=4;                   % cache α²
barrier_params.mu=50;% 2
barrier_params.p=r_obs;
barrier_params.k=100;


dummy_params = struct();
%DxArray = [100,500,1000,1500,2000];
DxArray = [10,200,500];
NR = 1;
filters = {
    SIRFilter(dummy_params, DxArray, NR),
    ENKFFilter(dummy_params, DxArray, NR),
    ENKF_Barrier_Filter(dummy_params, DxArray, NR),	
    APFFilter(dummy_params, DxArray, NR),
    BarrierFilter(dummy_params, DxArray, NR)


};

if is_N_fixed
  NArray=750*ones(size(DxArray));
else
   NArray=DxArray;
end







for iDX=1:numel(DxArray)
    Dx=DxArray(iDX);
    N=NArray(iDX);
    parhyper.r0=parhyper.Lm^(1/Dx);
   Dz = fix(3*Dx/5);          % set full observation
    %Dz=Dx;
    

    for nr = 1:NR 
        fixed_observed_components = randsample(Dx,Dz);
        fixed_observed_components= sort(fixed_observed_components);
        ok = 0;
        while not(ok)
    
            n_steps = ceil(5/he);%-----------------------------------------> up to half of final time ?
            Wx0 = sqrt(he)*randn([Dx n_steps]);% -------------------------> Brownian increment
            [x_ini,~] = exp_euler(rand([Dx 1]),he,F,n_steps,Dx,Wx0,sx);
            idx = randsample( fix(n_steps/2):n_steps, 1 );%------------------> choose an integer between nstep/2 ,nsptep
            x0 = x_ini(:,idx);
            % now we run with the 'regular' initialisation
            Wx = sqrt(he)*randn([Dx NTe]);%----------------------------------> create new increment
            [x,ok] = exp_euler(x0,he,F,NTe,Dx,Wx,sx); %--------------------> if ok =1, x does not  contain nan values?
                                                  % ------------------> x holds 10 000 realizations
                                                                      % of 800 state  
        end %while
        H0= eye(Dx) + randn([Dx Dx]).*5e-4;
        H0x = H0*x(1:Dx,(n_obs+1):n_obs:NTe+1); % ---------------->now skip the inital time,start from t=h ,take final time 
        ze_full = H0x + sz*randn(size(H0x));
        H=H0(fixed_observed_components,:); %%% CHOOSE SUB MATRİX AS THE OBSERVATION MATRIX
        ze_sparse = ze_full(fixed_observed_components,:); %------> observed values from t=h up to t=T
   

        X0 = x0 + sx*randn([Dx N]);
        %%% Ensure initial particles lie inside initial hypercube
        

        truth = x(:, filtered_solution_indices);
        truth_full= x;
        Pd_f = mean(sum(x(:, filtered_solution_indices).^2));
        Pd_p = mean(sum(x.^2));


        params = struct( ...
            'F',                      F, ...
            'sx',                     sx, ...
            'sz',                     sz, ...
            'he',                     he, ...
            'NTe',                    NTe, ...
            'n_obs',                  n_obs, ...
            'ze_sparse',              ze_sparse, ...
            'H',                      H, ...
            'X0',                     X0, ...
            'ness_thr',               ness_thr, ...
            'Dz',                     Dz, ...
            'fixed_observed_components', fixed_observed_components, ...
            'truth',                  truth, ...
            'truth_full',             truth_full, ...
            'Pd_f',                   Pd_f, ...
            'Pd_p',                   Pd_p, ...
            'parhyper',               parhyper, ...
            'barrier_params',         barrier_params );
        
       
        % Assign updated params to filters
        for k = 1:numel(filters)
            filters{k}.params = params;
            filters{k}.run(iDX, nr);
        end % end for filter
    end% end for number of experiment
end % end for Dx array


for k = 1:numel(filters)
    disp(filters{k}.metric);
end






f = figure('Visible', 'off');
hold on;

for k = 1:numel(filters)
    metric = filters{k}.metric;
    plot(metric.DxArray, mean(metric.MSEf, 2), '-o', 'DisplayName', metric.name);
end

xlabel('Dx');
ylabel('Mean MSE_f');
title('Mean Filtered MSE vs Dx');
legend('Location', 'best');
grid on;
savefig(f, 'mean_filtered_MSE_vs_Dx.fig');
close(f);
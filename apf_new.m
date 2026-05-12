function [Xf, Xp] = apf_new(F,sx,sz,h,NT,n_obs,z,H,X0,save_obs_indices,snapshot_cfg)
% function [Xf, Xp] = apf(F,sx,sz,h,NT,n_obs,z,H,X0)
% a
% F : Lorenz 96 forcing parameter
% sx : process noise std
% sz : observation noise std
% h : Euler integration step
% NT : no. of discrete time steps
% n_obs : observations are collected every n_obs discrete time units
% z : observations
% H : observation matrix
% X0 : initial particles 
% 
% Xf : filtered states
% Xp : predicted states
%
% Standard particle filter, resampling conditional on ESS
%

% recovers the no. of particles (N), the no. of slow oscillators (nosc)
[Dx, N] = size(X0);
if nargin < 10
    save_obs_indices = [];
end
if nargin < 11
    snapshot_cfg = struct();
end
save_obs_indices = unique(save_obs_indices(:))';

% initialisation
nt=NT/n_obs; % choose NT and obs such that obs| NT
Xf = zeros([Dx nt+1]);      % filtered estimates, at obs
Xp = zeros([Dx NT+1]);      % predicted estimates, at all times  
Xf(:,1) = mean(X0,2);
Xp(:,1) = mean(X0,2);
Xold = X0;                % auxiliary
w= ones([1 N])/N;        % initial weights

% noise variance
s2z = sz^2;
s2x=sx^2;
% time steps
sxsqrth=sx*sqrt(h);
save_snapshot_on_the_fly(snapshot_cfg, save_step_indices, 0, X0, w);

%
% --Time loop
%
for obs_idx=1:nt % it was i
     obs_z=z(:,obs_idx); %  i eliminated zeros, z(1) gives observed value at t=h 
    % auxiliary states & prediction
    Xa = Xold;
    Xr = Xold;
    
     for    inner_idx=1:n_obs 
        % auxiliary states
        Xdrift = l96dxdt(Xa,F,Dx);
        Xa = Xa + h*Xdrift;

        % predictions
        dWx=sxsqrth*randn(Dx,N);
        Xdrift = l96dxdt(Xr,F,Dx);
        Xr = Xr + h*Xdrift + dWx;
        %Xp(:,(i-1)*n_obs+jj+1) = mean(Xr,2);
        
        Xp(:,(obs_idx-1)*n_obs+inner_idx+1)=Xr*w'; % it was mean;
        
            
    end %jj

    % auxiliary weights
    allk = -(1/(2*s2z)) .* sum(( obs_z - H*Xa ).^2);  % log-likelihood
    alw = allk - max(allk);
    awu = exp(alw);
    aw = awu ./ sum(awu);

    % auxiliary resampling
    X = Xold;
    idxa = randsample(1:N, N, true, aw);
    X(:, 1:N) = X(:, idxa);

    % sampling up to the next observation
    for jj = 1:n_obs
        WX = sxsqrth*randn(Dx,N);
        Xdrift = l96dxdt(X,F,Dx);
        X = X + h*Xdrift + WX;
    end %jj    

    % weights
    llk = -(1/(2*s2z)) .* sum(( obs_z - H*X ).^2);  % log-likelihood
    lw=llk-allk(idxa);
    lw = lw - max(lw);
    wu = exp(lw);
    w = wu ./ sum(wu);
   
    % estimates

    Xf(:,obs_idx+1) = X*w';  %%% 
   % save_snapshot_on_the_fly(snapshot_cfg, save_obs_indices, obs_idx, X, w);
% --- NEW: Save snapshot at intermediate steps ---
        current_step = (obs_idx - 1) * n_obs + inner_idx;
        save_snapshot_on_the_fly(snapshot_cfg, save_obs_indices, current_step, Xnew, w);
    % resampling
    idx = randsample(1:N, N, true, w);
    Xold = X(:, idx);
    
end % observation index i

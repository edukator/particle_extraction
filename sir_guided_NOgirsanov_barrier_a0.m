function [Xf, Xp, resampling_counter] = sir_guided_NOgirsanov_barrier_a0(F,sx,sz,h,NT,n_obs,z,H,X0,ness_thr,barrier_params)

% Guided particle filter with baseline OU auxiliary AND a_n = 0.
%
% True model (particles):  dX = b(X) dt + sx dW
% Auxiliary model (for g): dZ = (-Z) dt + sx dW   (i.e., A=-I, a(t)=0)
%
% Guidance:
%   G(Delta, x) = ∇_x log g(t,x) = phi * H' * R(Delta)^{-1} ( y - H*(phi*x) )
% where
%   phi = exp(-Delta),
%   C(Delta) = (sx^2/2) (1 - exp(-2 Delta)) I,
%   R(Delta) = H C H' + sz^2 I = cvar*(H H') + sz^2 I.
%
% Weights over one window [t_{n-1}, t_n]:
%   log w <- log w + log g(Delta0, X_{t_{n-1}}) + ∫ < r(X_t), G(t,X_t) > dt
% with
%   r(x) = b(x) - (-x) = b(x) + x.

[Dx, N] = size(X0);
nt = NT / n_obs;  % require n_obs | NT

% Initialisation
Xf = zeros(Dx, nt+1);      % filtered mean at obs times
Xp = zeros(Dx, NT+1);      % predicted mean at all times
Xf(:,1) = mean(X0,2);
Xp(:,1) = mean(X0,2);

Xold = X0;
lw   = zeros(1, N);        % log weights
w    = ones(1, N)/N;       % weights

% Precompute
sxsqrth = sx*sqrt(h);
s2z=sz^2;

resampling_counter = 0;

HT      = H.';                       % transpose once

mu=barrier_params.mu;
p=barrier_params.p; % here is r0;
k=barrier_params.k;

prev_center=H*mean(X0,2);  % Initialize center of 1st hypercube in OBS from initial particle positions.

% -- Loop over observation windows
for obs_idx = 1:nt

    obs_z = z(:,obs_idx);
    C(:,1) = prev_center; %
    C(:,2) = obs_z;   %as Joaquin suggested;
    [ci,~]=get_centers_of_hypertube(h,n_obs,C);

    % ============================================================
    % (A) a_n = 0 (fixed)
    % ============================================================
    % No need to compute xbar, bbar, etc.

    % (B) log g at window start (Delta0 = n_obs*h)
    Delta0 = n_obs * h;
    logg0  = log_g_aux_OU_a0(Xold, obs_z, H, sx, sz, Delta0);  % 1 x N

    % accumulate integral: sum_k h <r(X_k), G_k>
    integ = zeros(1, N);

    % ============================================================
    % (C) Inner propagation with guided proposal
    % ============================================================
    for inner_idx = 1:n_obs

        % remaining time-to-observation
        Delta = (n_obs - inner_idx + 1) * h;

        % True drift b(Xold)
        Xdrift = l96dxdt(Xold, F, Dx);     % Dx x N (your vectorized function)


          % for barrier
        % barrier term is implemented 
        e    = ci(:,inner_idx) - H * Xold;                % d_y × N
      
        J = sum(e.^2, 1) ./ (2*s2z);  % 1xN  
       
        zeta = 1 ./ (1 + exp(-k * (J - p)));
         
        q     = HT * e; 
        gradL = -(1/s2z) * q .* zeta; 

        % Guidance G(Delta, Xold) = ∇ log g(t, Xold) under a_n = 0 auxiliary
        G = score_G_aux_OU_a0(Xold, obs_z, H, sx, sz, Delta);  % Dx x N

        % Remainder mismatch r(x) = b(x) - (-x) = b(x) + x
        rX = Xdrift + Xold;

        % Accumulate integral term
        integ = integ + h * sum(rX .* G, 1);

        % Guided Euler-Maruyama proposal:
        dWx  = sxsqrth * randn(Dx, N);

        Xdrift=Xdrift-mu*(sx^2)*gradL;% if you include barrier term to the weight computation
                                  % do this earlier
        
        Xnew = Xold + h * (Xdrift + (sx^2) * G) + dWx;

        % move forward
        Xold = Xnew;

        % predicted mean (store at every step)
        Xp(:, (obs_idx-1)*n_obs + inner_idx + 1) = Xnew * w.';
    end

    % ============================================================
    % (D) Weight update at observation time (guided)
    % ============================================================
   llk = -(1/(2*s2z)) .* sum(( obs_z - H*Xnew ).^2);  % log-likelihood
   lw = lw+ llk - max(lw+llk);
    wu = exp(lw);
    w  = wu ./ sum(wu);

    % filtered estimate at obs time
    Xf(:, obs_idx+1) = Xnew * w.';

    % ============================================================
    % (E) Resampling conditional on ESS
    % ============================================================
    NESS = (1/sum(w.^2))/N;
    if NESS < ness_thr
        idx = randsample(1:N, N, true, w);
        Xnew = Xnew(:, idx);
        Xold = Xnew;

        w  = ones(1, N)/N;
        lw = zeros(1, N);

        resampling_counter = resampling_counter + 1;
    end
        prev_center=C(:,2);

end
end

% ======================================================================
% log g under OU auxiliary with a_n = 0:
% Z_T | Z_t=x ~ N( phi*x, C(Δ) ), C(Δ)=cvar*I
% y = H Z_T + noise, noise~N(0, sz^2 I)
%
% log g = -0.5 * (y - H*(phi*x))^T R^{-1} (y - H*(phi*x)) + const
function logg = log_g_aux_OU_a0(X, y, H, sx, sz, Delta)
phi  = exp(-Delta);
cvar = (sx^2/2) * (1 - exp(-2*Delta));            % scalar so C=cvar*I
R    = cvar*(H*H.') + (sz^2)*eye(size(H,1));      % dy x dy

m = phi * X;                                      % Dx x N
innov = y - H*m;                                  % dy x N

U = R \ innov;                                    % dy x N
logg = -0.5 * sum(innov .* U, 1);                 % 1 x N
end

% Score G = ∇_x log g = phi * H' * R^{-1} (y - H*(phi*x))
function G = score_G_aux_OU_a0(X, y, H, sx, sz, Delta)
phi  = exp(-Delta);
cvar = (sx^2/2) * (1 - exp(-2*Delta));
R    = cvar*(H*H.') + (sz^2)*eye(size(H,1));

m = phi * X;                                      % Dx x N
innov = y - H*m;                                  % dy x N
U = R \ innov;                                    % dy x N

G = phi * (H.' * U);                              % Dx x N
end




















   
  

   

    
      

function [Xf, Xp] = sir_tamed(F,sx,sz,h,NT,n_obs,z,H,X0,ness_thr,save_obs_indices,snapshot_cfg,tamed_type)
% function [Xf, Xp] = sir(F,sx,sz,h,NT,n_obs,z,X0,ness_thr)
%
% F : Lorenz 96 forcing parameter
% sx : process noise std
% sz : observation noise std
% h : Euler integration step
% NT : no. of discrete time steps
% n_obs : observations are collected every n_obs discrete time units
% z : observations
% H : observation matrix
% X0 : initial particles 
% ness_thr : resampling threshold
% 
% Xf : filtered states
% Xp : predicted states
%
% Standard particle filter, resampling conditional on ESS
%

% recovers the no. of particles (N), the no. of slow oscillators (nosc)
[Dx, N] = size(X0);
if nargin < 11
    save_obs_indices = [];
end
if nargin < 12
    snapshot_cfg = struct();
end
save_obs_indices = unique(save_obs_indices(:))';

nt=NT/n_obs; % choose NT and obs such that obs| NT
% initialisation
Xf = zeros([Dx nt+1]);      % filtered estimates, at obs
Xp = zeros([Dx NT+1]);      % predicted estimates, at all times  
Xf(:,1) = mean(X0,2);
Xp(:,1) = mean(X0,2);
Xold = X0;                % auxiliary
lw = zeros([1 N]);        % log-weights  
w = ones([1 N])/N;        % weights  
% noise variance
s2z = sz^2;
% time steps
sxsqrth=sx*sqrt(h);

save_snapshot_on_the_fly(snapshot_cfg, save_obs_indices, 0, X0, w);
 
%
% --Time loop
%
%finer_idx_counter=1;
for obs_idx=1:nt

   obs_z=z(:,obs_idx); %  i eliminated zeros, z(1) gives observed value at t=h 
  

   %fprintf("observation index  of   %d  at time  , \n ",obs_idx);
   for inner_idx=1:n_obs
        
        Xdrift = l96dxdt(Xold,F,Dx);
        dWx=sxsqrth*randn(Dx,N);
        if strcmpi(tamed_type, 'Arnulf')
            Xdrift = Xdrift./(1+h*vecnorm(Xdrift)); % example taming function
            Xnew=  Xold + h*Xdrift+dWx;
        elseif strcmpi(tamed_type, 'ExpTamed')
            Xdrift = Xdrift.*  (1-exp(-h*abs(Xdrift) )  )./(abs(Xdrift)*h); % example taming function
            Xnew=  Xold + h*Xdrift+exp(-h*abs(Xdrift)).*dWx;
        else
            error('sir_tamed:UnknownTamedType', ...
                'Unknown tamed type "%s". Expected "Arnulf" or "ExpTamed".', ...
                string(tamed_type));
        end
        %finer_idx_counter=finer_idx_counter+1;

        % prediction
        Xold=Xnew;
        current_step = (obs_idx-1)*n_obs + inner_idx;
        Xp(:,current_step+1)=Xnew*w';
        save_snapshot_on_the_fly(snapshot_cfg, save_obs_indices, current_step, Xnew, w);
       % fprintf("finer_idx  % d, Xp stored at   %d \n", finer_idx_counter,(obs_idx-1)*n_obs+inner_idx+1);
   end 
   
   

   % prediction
     
       % Available observations

       % weights
       llk = -(1/(2*s2z)) .* sum(( obs_z - H*Xnew ).^2);  % log-likelihood
       lw = lw + llk - max(lw+llk);
       wu = exp(lw);
       w = wu ./ sum(wu);
       
       % estimation
       Xf(:,obs_idx+1) = Xnew*w';  %%% is it the correct place ?  (obs_idx+1)
       obs_step = obs_idx*n_obs;
       save_snapshot_on_the_fly(snapshot_cfg, save_obs_indices, obs_step, Xnew, w);
       
       % fprintf("finer_idx  % d, Xp stored at   %d \n", finer_idx_counter,(obs_idx-1)*n_obs+inner_idx+1);     
       % Resampling
       NESS = (1/sum(w.^2))/N; 
       if NESS<ness_thr
            idx = randsample(1:N, N, true, w);
            Xnew(:,1:N) = Xnew(:,idx);
           
            w = ones([1 N])/N;
            lw = zeros([1 N]);
       end %if
       
   %fprintf("--------------\n");
    
   % ...for the next time step
   Xold = Xnew;
    
end % time (n)

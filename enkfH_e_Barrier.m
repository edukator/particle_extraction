function [xf, xp] = enkfH_e_Barrier(F,sx,sz,h,NT,n_obs,z,H,X0,Dobs,obs_components,barrier_params)

% recovers the no. of particles (N),
[Dx, N] = size(X0);

% initialisation
Xold = X0;

nt=NT/n_obs; % choose NT and obs such that obs| NT

xf = zeros([Dx nt+1]);            % filtered estimates, 
xp = zeros([Dx n_obs*nt+1]);      % predicted estimates 


 
xp(:,1) = mean(X0,2);
xf(:,1) = mean(X0,2);

s2z = sz^2;
s2x=sx^2;
obs_counter=1; % to keep track of observation stored at only obs time, 

%% Barrier hyperparameters.
HT      = H.';                       % transpose once

mu=barrier_params.mu;
p=barrier_params.p; % here is r0;
k=barrier_params.k;

prev_center=H*mean(X0,2);  
C(:,1) = prev_center; %
C(:,2) = z(:,obs_counter);   %as Joaquin suggested;
[ci,ti]=get_centers_of_hypertube(h,n_obs,C); 
inner_idx=1;

for n = 2:(NT+1)  %   NOT  NT, IT IS   NT+1;
   

    if rem(n,n_obs) == 1
          % yeri burası değil.
         

          WX = sx*sqrt(h)*randn([Dx N]);
          Xdrift = l96dxdt(Xold,F,Dx);
          e    = ci(:,inner_idx) - H * Xold; 
          inner_idx=1;
          J = sum(e.^2, 1) ./ (2*s2z);
          zeta = 1 ./ (1 + exp(-k * (J - p)));
          q     = HT * e; 
          gradL = -(1/s2z) * q .* zeta; 
          Xdrift=Xdrift-mu*s2x*gradL;
          Xp = Xold + h*Xdrift + WX;

        % predictive mean
          xp(1:Dx,n) = mean(Xp,2);
           % predicted measurements & measurement covariance
          %Zp = Xp(obs_components,1:N);
          Zp=H*Xp;
          zp = mean(Zp,2);
          cZp = Zp-zp;
          Cz = cZp*cZp'./(N-1) + s2z*eye(Dobs);
              % cross-covariance
          cXp = Xp - xp(:,n);
          Cxz = (cXp*cZp')/(N-1);



           % Kalman gain
            try
                Kg = Cxz/Cz;
            catch ME
                fprintf(1,'EnKF w/ Euler: %s\n', ME.message);
                xf(:,n:NT) = nan;
                xp(:,n:NT) = nan;
                return;
            end% try end
              % updated ensemble
              %Zn = z(obs_idx,n) + sqrt(s2z)*randn([Dobs N]);
              
           Zn = z(:,obs_counter) + sqrt(s2z)*randn([Dobs N]);
           Xf = Xp + Kg*(Zn-Zp);

           % updated mean
           xf(:,obs_counter+1) = mean(Xf,2);


            % for the next step 
           Xold = Xf;
          prev_center=C(:,2); 
          obs_counter=obs_counter+1; 
          if n<NT+1 
            C(:,1) = prev_center; %
            C(:,2) = z(:,obs_counter);   %as Joaquin suggested;
            [ci,ti]=get_centers_of_hypertube(h,n_obs,C); 
          end 
          
          
          


          
    

     % end 
    else % if no observation
             % no observations: forward Euler sampling     

          WX = sx*sqrt(h)*randn([Dx N]);
          Xdrift = l96dxdt(Xold,F,Dx);
          e    = ci(:,inner_idx) - H * Xold; 
          inner_idx=inner_idx+1;
          J = sum(e.^2, 1) ./ (2*s2z);
          zeta = 1 ./ (1 + exp(-k * (J - p)));
          q     = HT * e; 
          gradL = -(1/s2z) * q .* zeta; 
           Xdrift=Xdrift-mu*s2x*gradL;
          %Xdrift=Xdrift./(1+h*vecnorm(Xdrift));
          Xp = Xold + h*Xdrift + WX;
          Xold = Xp;
          xp(1:Dx,n) = mean(Xp,2);% no need may be
          
           


    end % end if-else

   
      
end % end for


%fprintf("latest obs counter is %d",obs_counter-1);

end
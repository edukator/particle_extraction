
to run the code : main_cli(sx, sz, is_N_fixed, fig_dir)

for example;
main_cli(1.0, 1.0, true, 'results_run1')
----------------------------------------------
filename format is:

<filter_name>_obs_XXXX.mat (e.g., SIR_obs_0007.mat).f
-----------------------------------------
 inside snapshot, fields are:

obs_idx
Observation index where this snapshot was taken (e.g., 1, 7, 10).

particles
Particle/ensemble matrix at that observation instant.
Shape is typically Dx x N (state dimension × particle count).

weights
Weight vector for those particles (typically 1 x N).
For EnKF-style methods, this may be uniform weights. 

Dx 
State dimension 

------------------------------------------------------------------
To read the mat file 
S = load('snapshots/Dx_50/SIR_obs_0007.mat');
fieldnames(S.snapshot)
size(S.snapshot.particles)
size(S.snapshot.weights)
S.snapshot.obs_idx
S.snapshot.Dx
-------------------------------
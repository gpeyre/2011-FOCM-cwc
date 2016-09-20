% tests for various parameters of cwc : 
% smoothing of steps, hard/soft threshold, niter MCA

n = 512*2;
name = 'bv-steps';
rand('seed', 123456); randn('seed', 123456);

% params of medium
roughness = 3;
contrast = 3;
speed_blur = 5;
nrandomization = 10; % number of randomized frequencies
nsub = 10; % number of tested sub-sampling
nitermca = 1500*2;
save_iter = 1;
save_eps = 1;
additional_plot = 1;
add_str = '';
thresh = 'hard';
options.use_mad = 0;

sublist = linspace(.98,.02,nsub);

% use douglas rachford for L1 minimization
thresh = 'soft';
options.solver = 'dr';

test_wave_1d;
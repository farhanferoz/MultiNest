% example: estimate amplitude and initial phase of a
% sinusoidal signal in additive white gaussian noise

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

% SIMULATE DATA
% discrete times
dt = 1/100; % sampling rate of signal (Hz)
tlen = 3; % length of signal (seconds)
t = [0:dt:tlen]';

% sinusoidal signal model
amp = 10; % signal amplitude
phi = 2.3; % initial phase of signal (rads)
f = 24.788634; % signal frequency
s = sinusoid_model(t, {'amp', 'phi', 'f'}, {amp, phi, f});

% gaussian noise model (white)
sigma2 = 2; % variance of Gaussian noise
n = sqrt(sigma2) * randn(length(t), 1);

% data (t, y, covariance matrix)
y = s+n;
data{1} = t;
data{2} = y;
data{3} = sigma2;

clear amp phi s n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define nested sampling parameters
Nlive = 500;

Nmcmc = input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_gaussian;
model = @sinusoid_model;
prior = {'amp', 'uniform', 0, 20, 'fixed'; ...
         'phi', 'uniform', 0, 2*pi, 'cyclic'};
extraparams = {'f', f}; % fixed signal frequency

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc, 'diffevfrac', 2, 'covfrac', 10, 'walkfrac', 0, 'stretchfrac', 0);

% plot posterior distributions
wp = [1];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [2];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [1 2];
posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});


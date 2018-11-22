% example: a triangular likelihood function

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define nested sampling parameters
Nlive = 500;

Nmcmc = input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_model_likelihood;
model = @triangle_model;
prior = {'x', 'uniform', 0, 1, 'fixed'};
extraparams = {'gradient', 1}; % fixed signal frequency

data = 1; % dummy data

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);


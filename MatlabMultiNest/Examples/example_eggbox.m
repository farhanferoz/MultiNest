% example: an eggbox function

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
model = @eggbox_model;
prior = {'1', 'uniform', 0, 1, 'fixed'; ...
         '2', 'uniform', 0, 1, 'fixed'; ...
         '3', 'uniform', 0, 1, 'fixed'; ...
         '4', 'uniform', 0, 1, 'fixed'; ...
         '5', 'uniform', 0, 1, 'fixed'};
extraparams = {};

data = 1; % dummy data

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

% plot two of the dimensions and the likelihood
figure;
scatter3(post_samples(:,1), post_samples(:,2), post_samples(:,6), 20, ...
    post_samples(:,6), 'filled')

% example: (from hogg et al., 1008.4686)
% fit a line to data keeping outlier points

global verbose;
verbose = 1;
global DEBUG;
DEBUG = 0;

% read sample data for fitting a line
% NOTE: data{1} = x_i, data{2} = y_i, data{3} = sigma_yi
data = readdata_line;

% convert sigmas to a covariance matrix and reassign to data{3}
C = diag(data{3}.^2);
data{3} = C;

% define nested sampling parameters
Nlive = 500;
Nmcmc = input('Input number of iterations for MCMC sampling: (enter 0 for multinest sampling)\n');
tolerance = 0.1;
likelihood = @logL_gaussian;
model = @line_model;
prior = {'m', 'uniform', 0, 10, 'fixed'; ...
         'b', 'uniform', 0, 300, 'fixed'};
extraparams = {}; % no extra parameters beyond what's in the prior

% called nested sampling routine
[logZ, nest_samples, post_samples] = nested_sampler(data, Nlive, ...
  tolerance, likelihood, model, prior, extraparams, 'Nmcmc', Nmcmc);

% plot posterior distributions
wp = [1];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [2];
posteriors(post_samples, wp, {prior{wp,1}});
wp = [1 2];
posteriors(post_samples, wp, {prior{wp(1),1}, prior{wp(2),1}});


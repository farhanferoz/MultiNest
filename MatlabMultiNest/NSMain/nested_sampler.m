function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
          Nlive, tolerance, likelihood, model, prior, extraparams, ...
          varargin)

% function [logZ, nest_samples, post_samples] = nested_sampler(data, ...
%           Nlive, Nmcmc, tolerance, likelihood, model, prior, extraparams)
%
% This function performs nested sampling of the likelihood function from
% the given prior (given a set of data, a model, and a set of extra model
% parameters).
%
% By default the algorithm will draw new samples from a set of bounding
% ellipsoids constructed using the MultiNest algorithm for partitioning
% live points. However, if the optional 'Nmcmc' argument is set and
% Nmcmc > 0, new samples will be drawn from a proposal using an MCMC. This
% method is based on that of Veitch & Vecchio. For both methods the
% sampling will stop once the tolerance critereon has been reached.
%
% The likelihood should be the function handle of a likelihood function to
% use. This should return the log likelihood of the model parameters given
% the data.
%
% The model should be the function handle of the model function to be
% passed to the likelihood function.
%
% The prior should be a cell array with each cell containing five values:
%   parameter name (string)
%   prior type (string) e.g. 'uniform', 'gaussian' of 'jeffreys'
%   minimum value (for uniform prior), or mean value (for Gaussian prior)
%   maximum value (for uniform prior), or width (for Gaussian prior)
%   parameter behaviour (string):
%       'reflect' - if the parameters reflect off the boundaries
%       'cyclic'  - if the parameter space is cyclic
%       'fixed'   - if the parameters have fixe boundaries
%       ''        - for gaussian priors
%   e.g., prior = {'h0', 'uniform', 0, 1, 'reflect';
%                  'r', 'gaussian', 0, 5, '';
%                  'phi', 'uniform', 0, 2*pi, 'cyclic'};
%
% extraparams is a cell array of fixed extra parameters (in addition
% to those specified by prior) used by the model
% e.g.  extraparams = {'phi', 2;
%                      'x', 4};
%
% Optional arguments:
%  Set these via e.g. 'Nmcmc', 100
%   Nmcmc - if this is set then MultiNest will not be used as the sampling
%           algorithm. Instead an MCMC chain with this number of iterations
%           will be used to draw the number nested sample point.
%   Nsloppy - if this is set then during the MCMC the likelihood will only
%             be evaluted once every Nsloppy points rather than at every
%             iteration of the chain.
%   covfrac - the relative fraction of the iterations for which the MCMC
%             proposal distribution will be based on a Students-t
%             distribution defined by the covariance of the current live
%             points.
%   diffevfrac - the relative fraction of the iterations that will use
%                differential evolution to draw the new sample.
%   stretchfrac - the relative fraction of the iterations that will use the
%                 affine invariant ensemble stretch method for drawing a
%                 new sample
%   walkfrac - the relative fraction of the iterations that will use the
%              affine invariant ensemble walk method for drawing a new
%              sample
%   propscale - the scaling factor for the covariance matrix used by the
%               'covfrac' Students-t distribution proposal. This defaults
%               to 0.1.
%
% E.g. if covfrac = 10 then diffevfrac = 5 the Students-t proposal will be
% used 2/3s of the time and differential evolution 1/3. The default is to
% use the affine invariant samplers with the stretch move 75% of the time
% and the walk move 25% of the time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global verbose;
global DEBUG;

Nmcmc = 0; % default with this set to zero is to use MultiNest
Nsloppy = 0;
covfrac = 0;
diffevfrac = 0;
walkfrac = 25;
stretchfrac = 75;
propscale = 0.1;

% get optional input arguments
optargin = size(varargin,2);

% get optional arguments
if optargin > 1
    for i = 1:2:optargin
        if strcmpi(varargin{i}, 'Nmcmc') % number of MCMC samples
            if ~isempty(varargin{i+1})
                if varargin{i+1} < 1
                    fprintf(1, 'Using MultiNest algorithm\n');
                else
                    Nmcmc = varargin{i+1};
                end
            end
        elseif strcmpi(varargin{i}, 'Nsloppy') % number of burn in samples
            if ~isempty(varargin{i+1})
                if varargin{i+1} < 0
                    fprintf(1, 'Number of \"sloppy\" samples is silly. Setting to zero\n');
                else
                    Nsloppy = varargin{i+1};
                end
            end
        elseif strcmpi(varargin{i}, 'covfrac') % fraction of MCMC iterations using Students't proposal
            if varargin{i+1} > 0
                covfrac = varargin{i+1};
            end
        elseif strcmpi(varargin{i}, 'diffevfrac') % fraction of MCMC iterations using differential evolution
            if varargin{i+1} > 0
                diffevfrac = varargin{i+1};
            end
        elseif strcmpi(varargin{i}, 'walkfrac') % fraction of MCMC iterations using walk move
            if varargin{i+1} >= 0
                walkfrac = varargin{i+1};
            end
        elseif strcmpi(varargin{i}, 'stretchfrac') % fraction of MCMC iterations using stretch move
            if varargin{i+1} >= 0
                stretchfrac = varargin{i+1};
            end
        elseif strcmpi(varargin{i}, 'propscale') % the scaling factor for the covariance matrix
            if varargin{i+1} > 0
               propscale = varargin{i+1};
            end
        end
    end
end

% get the number of parameters from the prior array
D = size(prior,1);

% get all parameter names
parnames = prior(:,1);

if ~isempty(extraparams)
    extraparnames = extraparams(:,1);
    extraparvals = extraparams(:,2);
    parnames = cat(1, parnames, extraparnames);
else
    extraparvals = [];
end

% draw the set of initial live points from the prior
livepoints = zeros(Nlive, D);

for i=1:D
    priortype = prior{i,2};
    p3 = prior{i,3};
    p4 = prior{i,4};

    % currently only handles uniform or Gaussian priors
    if strcmp(priortype, 'uniform')
        livepoints(:,i) = p3 + (p4-p3)*rand(Nlive,1);
    elseif strcmp(priortype, 'gaussian')
        livepoints(:,i) = p3 + p4*randn(Nlive,1);
    elseif strcmp(priortype, 'jeffreys')
        % uniform in log space
        livepoints(:,i) = 10.^(log10(p3) + (log10(p4)-log10(p3))*rand(Nlive,1));
    end
end

% check whether likelihood is a function handle, or a string that is a
% function name
if ischar(likelihood)
    flike = str2func(likelihood);
elseif isa(likelihood, 'function_handle')
    flike = likelihood;
else
    error('Error... Expecting a model function!');
end

% calculate the log likelihood of all the live points
logL = zeros(Nlive,1);

for i=1:Nlive
    parvals = cat(1, num2cell(livepoints(i,:)'), extraparvals);
    logL(i) = feval(flike, data, model, parnames, parvals);
end

% now scale the parameters, so that uniform parameters range from 0->1,
% and Gaussian parameters have a mean of zero and unit standard deviation
for i=1:Nlive
    livepoints(i,:) = scale_parameters(prior, livepoints(i,:));
end

% initial tolerance
tol = inf;

% initial width of prior volume (from X_0=1 to X_1=exp(-1/N))
logw = log(1 - exp(-1/Nlive));

% initial log evidence (Z=0)
logZ = -inf;

% initial information
H = 0;

% initialize array of samples for posterior
nest_samples = zeros(1,D+1);

%%%%%%%%%%%%%%%%
% some initial values if MultiNest sampling is used
h = 1.1; % h values from bottom of p. 1605 of Feroz and Hobson
FS = h; % start FS at h, so ellipsoidal partitioning is done first time
K = 1; % start with one cluster of live points

% get maximum likelihood
logLmax = max(logL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize iteration counter
j = 1;

%figure;

% MAIN LOOP
while tol > tolerance || j <= Nlive

    % expected value of true remaining prior volume X
    VS = exp(-j/Nlive);

    % find minimum of likelihoods
    [logLmin, idx] = min(logL);

    % set the sample to the minimum value
    nest_samples(j,:) = [livepoints(idx, :) logLmin];

    % get the log weight (Wt = L*w)
    logWt = logLmin + logw;

    % save old evidence and information
    logZold = logZ;
    Hold = H;

    % update evidence, information, and width
    logZ = logplus(logZ, logWt);
    H = exp(logWt - logZ)*logLmin + ...
        exp(logZold - logZ)*(Hold + logZold) - logZ;
    %logw = logw - logt(Nlive);
    logw = logw - 1/Nlive;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Nmcmc > 0

        % do MCMC nested sampling

        % get the Cholesky decomposed covariance of the live points
        % (do every 100th iteration - CAN CHANGE THIS IF REQUIRED)
        if mod(j-1, 100) == 0
            % NOTE that for numbers of parameters >~10 covariances are often
            % not positive definite and cholcov will have "problems".
            %cholmat = cholcov(propscale*cov(livepoints));

            % use modified Cholesky decomposition, which works even for
            % matrices that are not quite positive definite
            % from http://infohost.nmt.edu/~borchers/ldlt.html
            % (via http://stats.stackexchange.com/questions/6364
            % /making-square-root-of-covariance-matrix-positive-definite-matlab
            cv = cov(livepoints);
            [l, d] = mchol(propscale*cv);
            cholmat = l.'*sqrt(d);

            %plot3(livepoints(:,1), livepoints(:,2), livepoints(:,3), 'r.');
            %drawnow();
        end

        % draw a new sample using mcmc algorithm
        [livepoints(idx, :), logL(idx)] = draw_mcmc(livepoints, cholmat, ...
              logLmin, prior, data, flike, model, Nmcmc, Nsloppy, ...
              covfrac, diffevfrac, walkfrac, stretchfrac, parnames, ...
              extraparvals);

    else

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % do MultiNest nested sampling

        % separate out ellipsoids
        if FS >= h
            % NOTE: THIS CODE IS GUARANTEED TO RUN THE 1ST TIME THROUGH
            % calculate optimal ellipsoids
            [Bs, mus, VEs, ns] = optimal_ellipsoids(livepoints, VS);
            K = length(VEs); % number of ellipsoids (subclusters)

        else
            % simply rescale the bounding ellipsoids
            for k=1:K
                scalefac = max([1 (exp(-(j+1)/Nlive)*ns(k)/Nlive)/VEs(k)]);

                % scale bounding matrix and volume
                if scalefac ~= 1
                    Bs((k-1)*D+1:k*D,:) = Bs((k-1)*D+1:k*D,:)*scalefac^(2/D);
                    VEs(k) = scalefac*VEs(k);
                end
           end

        end

        if DEBUG && D==2
           % plot 2-dimensionsal live points and bounding ellipses
           plot_2d_livepoints_with_ellipses(livepoints, Bs, mus);
        end

        % calculate ratio of volumes (FS>=1) and cumulative fractional volume
        Vtot = sum(VEs);
        FS = Vtot/VS;
        fracvol = cumsum(VEs)/Vtot;

        % draw a new sample using multinest algorithm
        [livepoints(idx, :), logL(idx)] = draw_multinest(fracvol, ...
              Bs, mus, logLmin, prior, data, flike, model, ...
              parnames, extraparvals);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update maximum likelihood if appropriate
    if logL(idx) > logLmax
        logLmax = logL(idx);
    end

    % work out tolerance for stopping criterion
    tol = logplus(logZ, logLmax - (j/Nlive)) - logZ;

    % display progress (optional)
    if verbose
        fprintf(1, 'log(Z): %.5e, tol = %.5e, K = %d, iteration = %d\n', ...
                logZ, tol, K, j);
    end

    % update counter
    j = j+1;

end

% sort the remaining points (in order of likelihood) and add them on to
% the evidence
[logL_sorted, isort] = sort(logL);
livepoints_sorted = livepoints(isort, :);

for i=1:Nlive
    logZ = logplus(logZ, logL_sorted(i) + logw);
end

% append the additional livepoints to the nested samples
nest_samples = [nest_samples; livepoints_sorted logL_sorted];

% rescale the samples back to their true ranges
for i=1:length(nest_samples)
    nest_samples(i,1:end-1) = ...
     rescale_parameters(prior, nest_samples(i,1:end-1));
end

% convert nested samples into posterior samples - nest2pos assumes that the
% final column in the sample chain is the log likelihood
post_samples = nest2pos(nest_samples, Nlive);

return

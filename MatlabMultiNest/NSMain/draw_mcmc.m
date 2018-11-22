function [sample, logL] = draw_mcmc(livepoints, cholmat, logLmin, ...
    prior, data, likelihood, model, Nmcmc, Nsloppy, covfrac, ...
    diffevfrac, walkfrac, stretchfrac, parnames, extraparvals)

% function [sample, logL] = draw_mcmc(livepoints, cholmat, logLmin, ...
%    prior, data, likelihood, model, Nmcmc, Nsloppy, covfrac, ...
%    diffevfrac, walkfrac, stretchfrac, parnames, extraparvals)
%
% This function will draw a multi-dimensional sample from the prior volume
% for use in the nested sampling algorithm. The new point will have a
% likelihood greater than the value logLmin. The new point will be found by
% evolving a random multi-dimensional sample from within the sample array,
% livepoints, using an MCMC with Nmcmc iterations.
%
% The MCMC can use four different proposals to draw new samples:
%   - a Students-t (with N=2 degrees of freedon) proposal distribution
%     based on the Cholesky decomposed covariance matrix of the array,
%     cholmat.
%   - using differential evolution by taking two random points from the
%     current live points.
%   - using the affine invariant stretch move
%   - using the affine invariant walk move
%     (see Goodman & Weare (2010) DOI: 10.2140/camcos.2010.5.65)
% Each of these proposals will be used in propotion to the values given by
% covfrac, diffevfrac, stretchfrac and walkfrac respectively.
%
% extraparvals is a vector of additional parameters needed by the model.
%
% If Nsloppy is an integer greater than 1, then the likelihood will only be
% evaluated once every Nsloppy point, otherwise samples will just be
% accepted/rejected based on the prior.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global verbose;

l2p = 0.5*log(2*pi); % useful constant
maxscale = 3;
lms = log(maxscale);
walksamplesize = 3;
totalfrac = covfrac+diffevfrac+walkfrac+stretchfrac;

Nlive = size(livepoints,1);
Npars = size(livepoints,2);

if walkfrac > 0
    assert(walksamplesize <= Nlive, 'Error... number of livepoints must be greater than 3')
end

Ndegs = 2; % degrees of freedom of Students't distribution

% initialize counters
acctot = 0;
Ntimes = 1;

while 1
    acc = 0;

    % get random point from live point array
    sampidx = ceil(rand*Nlive);
    sample = livepoints(sampidx, :);

    % get the sample prior
    currentPrior = 0;

    for j=1:Npars
        priortype = prior{j,2};
        p3 = prior{j,3};
        p4 = prior{j,4};

        if strcmp(priortype, 'uniform')
            pv = -log(p4-p3);
            currentPrior = currentPrior + pv;
        elseif strcmp(priortype, 'gaussian')
            pv = -l2p - sample(j)^2/2;
            currentPrior = currentPrior + pv;
        elseif strcmp(priortype, 'jeffreys')
            pv = -log(10^(sample(j)*(log10(p4) - log10(p3)) + log10(p3)));
            currentPrior = currentPrior + pv;
        end
    end

    for i=1:Nmcmc
        propratio = 0; % proposal ratio
        rnum = rand*totalfrac;
        if rnum < covfrac % use Students-t proposal with covariance matrix
            % draw points from mulitvariate Gaussian distribution
            gasdevs = randn(Npars,1);
            sampletmp = (cholmat*gasdevs)';

            % calculate chi-square distributed value
            chi = sum(randn(Ndegs,1).^2);

            % add value onto old sample
            sampletmp = sample + sampletmp*sqrt(Ndegs/chi);
        elseif rnum < (covfrac+diffevfrac) % use differential evolution
            % use differential evolution
            % draw two random (different points) A and B and add (B-A) to
            % the current sample
            idx1 = ceil(rand*Nlive);
            idx2 = idx1;
            while idx2 == idx1
                idx2 = ceil(rand*Nlive);
            end

            A = livepoints(idx1, :);
            B = livepoints(idx2, :);

            sampletmp = sample + (B-A);
        elseif rnum < (covfrac+diffevfrac+stretchfrac) % use an ensemble sampler (with Walk and Stretch moves)
            % stretch move
            R = rand;
            logscale = 2.*lms*R - lms;
            scale = exp(logscale); % scale factor
            propratio = logscale*Npars; % proposal ratio

            % pick two random samples from the live points
            idx1 = randi(Nlive);
            idx2 = randi(Nlive);
            while idx1 == idx2
                idx2 = randi(Nlive);
            end

            % get new sample
            A = livepoints(idx1, :);
            B = livepoints(idx2, :);

            sampletmp = sample + scale*(B-A);
        else
            % otherwise use walk move
            idxs = randperm(Nlive, walksamplesize);

            % get centre of mass of points
            com = mean(livepoints(idxs,:));

            % generate univariate random Gaussian variables for each sample
            uv = randn(walksamplesize,1);

            step = (livepoints(idxs,:)-repmat(com,walksamplesize,1))'*uv;

            sampletmp = sample + step.';
        end

        % check sample is within the (scaled) prior
        newPrior = 0;
        for j=1:Npars
            priortype = prior{j,2};
            p3 = prior{j,3};
            p4 = prior{j,4};

            if strcmp(priortype, 'uniform')
                behaviour = prior{j,5};

                dp = 1;

                if sampletmp(j) < 0 || sampletmp(j) > 1
                    if strcmp(behaviour, 'reflect')
                        % reflect the value from the boundary
                        sampletmp(j) = 1 - mod(sampletmp(j), dp);
                    elseif strcmp(behaviour, 'cyclic')
                        % wrap parameter from one side to the other
                        while sampletmp(j) > 1
                            sampletmp(j) = sampletmp(j) - 1;
                        end

                        while sampletmp(j) < 0
                            sampletmp(j) = sampletmp(j) + 1;
                        end
                    else
                        newPrior = -inf;
                        break;
                    end

                end

                pv = -log(p4-p3);
                newPrior = newPrior + pv;

            elseif strcmp(priortype, 'gaussian')
                pv = -l2p - sampletmp(j)^2/2;
                newPrior = newPrior + pv;
            elseif strcmp(priortype, 'jeffreys')
                behaviour = char(prior(j,5));

                dp = 1;

                if sampletmp(j) < 0 || sampletmp(j) > 1
                    if strcmp(behaviour, 'reflect')
                        % reflect the value from the boundary
                        sampletmp(j) = 1 - mod(sampletmp(j), dp);
                    else
                        newPrior = -inf;
                        break;
                    end
                end

                pv = -log(10^(sampletmp(j)*(log10(p4) - log10(p3)) + log10(p3)));
                newPrior = newPrior + pv;
            end
        end

        if log(rand) > newPrior - currentPrior + propratio % reject point
            continue;
        elseif Nsloppy ~= 0 && mod(i, Nsloppy) ~= 0 % don't need to recalculate the posterior at every point
            sample = sampletmp;
            currentPrior = newPrior;
            continue;
        end

        % rescale sample back to its proper range for likelihood
        sc = rescale_parameters(prior, sampletmp);

        % get the likelihood of the new sample
        %likestart = tic;
        logLnew = feval(likelihood, data, model, parnames, ...
                        cat(1, num2cell(sc), extraparvals));
        %likedur = toc(likestart);
        %fprintf(1, 'liketime = %.6f\n', likedur);

        % if logLnew is greater than logLmin accept point
        if logLnew > logLmin
            acc = acc + 1;
            currentPrior = newPrior;
            sample = sampletmp;
            logL = logLnew;
        end
    end

    % only break if at least one point was accepted otherwise try again
    if acc > 0
        acctot = acc;
        break;
    end

    acctot = acctot + acc;
    Ntimes = Ntimes + 1;
end

% print out acceptance ratio
if verbose
    if ~Nsloppy
        fprintf(1, 'Acceptance ratio: %1.4f, ', acctot/(Ntimes*Nmcmc));
    else
        fprintf(1, 'Acceptance ratio: %1.4f, ', acctot*Nsloppy/(Ntimes*Nmcmc));
    end
end

return

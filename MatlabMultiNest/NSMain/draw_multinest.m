function [sample, logL] = draw_multinest(fracvol, Bs, mus, ...
    logLmin, prior, data, likelihood, model, parnames, extraparvals)

% function [sample, logL] = draw_multinest(fracvol, Bs, mus, ...
%     logLmin, prior, data, likelihood, model, parnames, extraparvals)
%
% This function draws a multi-dimensional sample from the prior volume
% for use in the nested sampling algorithm. The new point will have a
% likelihood greater than the value logLmin. The new point will be found by
% drawing a random multi-dimensional sample from within the set of optimal
% ellipsoids constructed using the MultiNest algorithm.  The bounding
% ellipsoids are defined by their bounding matrices Bs and centroids mus.
% extraparvals is a vector of additional parameters needed by the model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG

% extra number of ellipsoids, number of dimensions
K = size(mus, 1);
ndims = size(mus,2);

while 1
    % find the ellipsoid from which to draw a new point
    rval = rand;

    for k=1:K
        if rval < fracvol(k)
            continue
        else
            break
        end
    end
    k0 = k;

    % extract bounding matrix and centroid for that ellipsoid
    B = Bs((k0-1)*ndims+1:k0*ndims,:);
    mu = mus(k0,:);

    % draw points from that ellipsoid until logL >= logLmin
    logL = -inf;
    while logL < logLmin
        in_range = 1; % default value

        % draw one point from the ellipsoid
        pnt = draw_from_ellipsoid(B, mu, 1);

        % make sure that the point lies in unit hypercube
        for ii=1:ndims
             if pnt(ii)<0 || pnt(ii)>1
                 in_range = 0;
                 if DEBUG; fprintf('new point not in range!!!!\n'); end
             end
        end

        if in_range
            % assign as candidate replacement live point
            sample = pnt;

            % rescale point back to full range
            rescaledpnt = rescale_parameters(prior, pnt);

            % get new likelihood
            logL = feval(likelihood, data, model, parnames, ...
                         cat(1, num2cell(rescaledpnt), extraparvals));
        end

    end

    % check how many ellipsoids this point lies in
    inN = in_ellipsoids(pnt, Bs, mus);

    % only accept sample with 1/inN probability
    if rand < 1/inN
        break
    end
end

return


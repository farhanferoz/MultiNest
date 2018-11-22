function [Bs, mus, VEs, ns] = optimal_ellipsoids(u, VS)

% function [Bs, mus, VEs, ns] = optimal_ellipsoids(u, VS)
%
% This function attempts to optimally partition the multi-dimensional
% samples u (uniformly distributed within the sample volume VS), into
% a set of subclusters enclosed by bounding ellipsoids.  The algorithm
% is based on Algorithm 1 of the MULTINEST paper by Feroz, Hobson,
% and Bridges, MNRAS, 398, 1601-1614 (2009).
%
% Output:
%   Bs:  an array of bounding matrices for the ellipsoids enclosing
%        the subclusters, scaled to have at least the minimum volume
%        required by the subclusters. ( (K x ndims) x ndims )
%   mus: an array of centroids for the bounding ellipsoids (K x ndims)
%   VEs: an array of volumes for the bounding ellipsoids   (K x 1)
%   ns:  an array containing the number of points for each subcluster (K x 1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG

N = size(u,1); % number of samples in multi-dimensional space
ndims = size(u,2); % number of dimensions

% calculate bounding matrix, etc. for bounding ellipsoid associated
% with the original set of points u
[B, mu, VE, flag] = calc_ellipsoid(u, VS);

% attempt to split u into two subclusters
[u1, u2, VE1, VE2, nosplit] = split_ellipsoid(u, VS);
n1 = size(u1,1);
n2 = size(u2,1);

if nosplit || n1<ndims+1 || n2<ndims+1
    % couldn't split the cluster
    Bs = B;
    mus = mu;
    VEs = VE;
    ns = N;
else
    % check if we should keep the partitioning of S
    if (VE1 + VE2 < VE || VE > 2*VS)
        if DEBUG
            fprintf('PARTITION ACCEPTED: N=%d splits to n1=%d, n2=%d\n', N, n1, n2);
        end
        VS1 = n1 * VS / N;
        VS2 = n2 * VS / N;

        [B1, mu1, VE1, n1] = optimal_ellipsoids(u1, VS1);
        [B2, mu2, VE2, n2] = optimal_ellipsoids(u2, VS2);

        Bs =  [B1 ; B2];
        mus = [mu1 ; mu2];
        VEs = [VE1 ; VE2];
        ns =  [n1 ; n2];
    else
        if DEBUG
            fprintf('PARTITION REJECTED: N=%d doesnt split into n1=%d and n2=%d\n', N, n1, n2);
        end
        Bs = B;
        mus = mu;
        VEs = VE;
        ns = N;
    end

end

return


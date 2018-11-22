function N = in_ellipsoids(pnt, Bs, mus)

% function N = in_ellipsoids(pnt, Bs, mus)
%
% This function works out how many of the ellipsoids (defined by the
% bounding matrices Bs and centroids mus) contain the point pnt.
% This number is returned in N.
%
% Bs is a [(Kxndims) x ndims] array, where K=total number of ellipsoids
% and ndims = dimension of the parameter space.
% mus is a [K x ndims] array.
% pnt is a ndims-dimensional vector.
%
% NOTE: in the future it may be quicker to input precalculated eigenvalues
% and eigenvectors into this function rather than the bounding matrices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 0;

% total number of ellipsoids and number of dimensions
K = size(mus,1);
ndims = size(mus,2);

% loop over number of ellipsiods and work out whether it contains the point
for k=1:K
    % set the point to have the same origin as the ellipsoid
    pntnew = pnt - mus(k,:);

    % extract the bounding matrix
    B = Bs((k-1)*ndims+1:k*ndims,:);

    % get the eigenvalues and eigenvectors of the ellipsoid
    [V, E] = eig(B); % V is matrix of eigenvectors (as columns)
    D = sqrt(diag(E));

    % rotate points to be on coordinate axes of the ellipsiod
    pntnew = pntnew * V;

    % scale points so that it's equivalent to having unit hyper-spheroids
    % rather than ellipsiods
    pntnew = pntnew ./ D';

    % get distance to point from centre of hyper-sphere
    dpnt = sqrt(sum(pntnew.^2));

    if dpnt <= 1
        % values is within the ellipsiod
        N = N + 1;
    end
end

return


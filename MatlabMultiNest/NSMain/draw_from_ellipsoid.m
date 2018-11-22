function pnts = draw_from_ellipsoid(B, mu, N )

% function pnts = draw_from_ellipsoid(B, mu, N )
%
% This function draws points uniformly from an ndims-dimensional ellipsoid
% with edges and orientation defined by the the bounding matrix B and
% centroid mu.  The output is a Nxndims dimensional array pnts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get number of dimensions from the bounding matrix B
ndims = size(B,1);

% calculate eigenvalues and vectors of the bounding matrix
[V, E] = eig(B);
D = sqrt(diag(E));

% check size of mu and transpose if necessary
if size(mu,1) > 1
    mu = mu';
end

% generate radii of hyperspheres
rs = rand(N,1);

% generate points
pt = randn(N,ndims);

% get scalings for each point onto the surface of a unit hypersphere
fac = sum(pt.^2, 2);

% calculate scaling for each point to be within the unit hypersphere
% with radii rs
fac = (rs.^(1/ndims)) ./ sqrt(fac);

pnts = zeros(N,ndims);

% scale points to the ellipsoid using the eigenvalues and rotate with
% the eigenvectors and add centroid
for i=1:N
    % scale points to a uniform distribution within unit hypersphere
    pnts(i,:) = fac(i)*pt(i,:);

    % scale and rotate to ellipsoid
    pnts(i,:) = (pnts(i,:) .* D' * V') + mu;
end

return

function [B, mu, VE, flag] = calc_ellipsoid(u, VS)
%
% calculate properties of ellipsoid given a set of points u
%
% Inputs:
%    u:  Nxndims array where N is the number point and ndims is the
%        number of dimensions
%    VS: minimum volume that the bounding ellipsoid should have
%
% Outputs:
%    B:    bounding matrix for ellipsoid including scale factor
%          for mininimum volume
%    mu:   centroid
%    VE:   volume of ellipsoid
%    flag: = 1 if number of points too small or bounding matrix
%          has bad condition number; otherwise = 0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG

% default values
B = [];
mu = [];
VE = [];
flag = 0;

% extract number of points and number of dimensions
N = size(u, 1);
ndims = size(u, 2);

% check that total number of points is large enough
if N < ndims+1
    if DEBUG; fprintf('number of samples too small to calculate bounding matrix for ellipsoid\n'); end;
    flag = 1;
    return;
end

% constant factor for volume of ellipsoid
const = pi^(ndims/2)/gamma(ndims/2 + 1);

% calculate covariance matrix and centroid
C = cov(u);
mu = mean(u);

% check condition number of C (eps = 2.2204e-16)
if rcond(C)<eps || isnan(rcond(C))
    if DEBUG; fprintf('bad condition number!\n'); end;
    flag = 1;
    return;
end

% find scale factor for bounding ellipsoid E
fB = 0;
for i=1:N
    f = ( (u(i,:)-mu) / C ) * (u(i,:)-mu)';
    if f > fB
        fB = f;
    end
end

% calculate volume of bounding ellipsoid E
VE = const*sqrt(det( fB * C ));

% expand volume of bounding ellipsoid to VS if necessary
fV = 1;
if VE < VS
    fV = (VS/VE)^(2/ndims);
    VE = VS;
end

% scale C to get bounding matrix B
B = fV * fB * C;

return

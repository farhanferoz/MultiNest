function [u1, u2, VE1, VE2, nosplit] = split_ellipsoid(u, VS)

% function [u1, u2, VE1, VE2, nosplit] = split_ellipsiod(u, VS)
%
% This function takes in a set of multi-dimensional data points u and the
% sample volume (VS) that they occupy. It uses the k-means algorthim to
% split the points into two sub-clusters and uses an optimisation scheme to
% re-assign points, if necessary, between the sub-clusters. This is based
% on the description in Algorithm 1 of the MULTINEST paper by Feroz,
% Hobson, and Bridges, MNRAS, 398, 1601-1614 (2009).
%
% The function returns the points in the two sub-cluster u1 and u2, and
% the volumes of the ellipsoid subclusters VE1 and VE2.  The flag nosplit
% is set to 1 if the splitting cannot be done; otherwise = 0.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global DEBUG
max_attempt = 50; % maximum number of attempts to recluster points

% default return values
nosplit = 0;
VE1 = [];
VE2 = [];
u1 = [];
u2 = [];

% extract number of samples and number of dimensions
N = size(u,1);
D = size(u,2);

% check total number of samples
if N < 2*(D+1)
    if DEBUG; fprintf('CANT SPLIT: total number of samples is too small!  N = %d\n', N); end;
    nosplit = 1;
    return;
end

% use kmeans to separate the data points into two sub-clusters
[idx, mu] = kmeans(u,2);
u1 = u(idx==1,:);
u2 = u(idx==2,:);

n1 = size(u1,1); % number of samples in S1
n2 = size(u2,1); % number of samples in S2

% check number of points in subclusters
if n1 < D+1 || n2 < D+1
    if DEBUG; fprintf('CANT SPLIT: number of samples in subclusters is too small! n1 = %d, n2 = %d\n', n1, n2); end;
    nosplit = 1;
    return;
end

% preallocate temp arrays
temp_u1 = cell(max_attempt,1);
temp_u2 = cell(max_attempt,1);
temp_VE1 = zeros(max_attempt,1);
temp_VE2 = zeros(max_attempt,1);
FS = zeros(max_attempt,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
numreassigned = [];
counter = 1;
while 1

    % calculate minimum volume of ellipsoids
    VS1 = VS*n1/N;
    VS2 = VS*n2/N;

    % calculate properties of bounding ellipsoids for the two subclusters
    [B1, mu1, VE1, flag1] = calc_ellipsoid(u1, VS1);
    [B2, mu2, VE2, flag2] = calc_ellipsoid(u2, VS2);

    % check flags
    if flag1 || flag2
        if DEBUG; fprintf('CANT SPLIT!!\n'); end;
        nosplit = 1;
        break
    end

    % construct temporary arrays and cell arrays containing results for
    % each pass through the loop
    temp_u1{counter} = u1;
    temp_u2{counter} = u2;
    temp_VE1(counter) = VE1;
    temp_VE2(counter) = VE2;
    FS(counter)=(VE1+VE2)/VS;

    % DEBUG print statement
    if DEBUG
        fprintf('SPLIT ELLIPSOID: counter = %d, FS = %f, numreassigned = %d\n', ...
                counter, (VE1+VE2)/VS, numreassigned);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check if points need to be reassigned to the other subcluster
    reassign = 0;

    m1 = 0;
    m2 = 0;
    u1new = zeros(N,D);
    u2new = zeros(N,D);

    % for all points get the Mahalanobis distance between each point and
    % the centroid of each ellipse and assign accordingly
    numreassigned = 0;
    for i=1:N
        % get d = (u-mu)^T * B^-1 * (u-mu)
        du1 = ( (u(i,:)-mu1) / B1 ) * (u(i,:)-mu1)';
        du2 = ( (u(i,:)-mu2) / B2 ) * (u(i,:)-mu2)';

        % calculate hk = VEk * duk / VSk;
        h1 = VE1 * du1 / VS1;
        h2 = VE2 * du2 / VS2;

        if h1 < h2
            m1 = m1 + 1;
            u1new(m1, :) = u(i,:);

            % check if point has been reassigned or not
            if idx(i) ~= 1
                reassign = 1;
                idx(i) = 1;
                numreassigned = numreassigned +1;
            end
        else
            m2 = m2 + 1;
            u2new(m2, :) = u(i,:);

            % check if point has been reassigned or not
            if idx(i) ~= 2
                reassign = 1;
                idx(i) = 2;
                numreassigned = numreassigned +1;
            end

        end

    end

    clear u1 u2 mu1 mu2;
    n1 = m1;
    n2 = m2;

    u1 = u1new(1:n1,:);
    u2 = u2new(1:n2,:);

    clear u1new u2new;

    % update counter
    counter = counter + 1;

    if reassign && counter <= max_attempt
        continue
    else
        % DEBUG print statement
        if DEBUG
            fprintf('SPLIT ELLIPSOID: counter = %d, FS = %f, numreassigned = %d\n', counter, (VE1+VE2)/VS, numreassigned);
            if counter > max_attempt
                fprintf('SPLIT ELLIPSOID: exceeded maximum attempts; take min F(S).\n');
            end
        end

        break
    end

end

% find minimum F(S) and return
[minFS, idx] = min(FS);
u1 = temp_u1{idx};
u2 = temp_u2{idx};
VE1 = temp_VE1(idx);
VE2 = temp_VE2(idx);

if DEBUG; fprintf('SPLIT ELLIPSOID: min F(S) = %f\n', minFS); end;

return


% script to test optimal_ellipsoids.m program
close all
clear all
rand('state',1);
randn('state',1);
global DEBUG
%DEBUG = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 1: single ellipsoid, 2d
clear B mu Bs mus x y z data
D = 2;

% generate uniform samples from a single ellipsoid
mu = [0 0];
B = [3 .2; .2 0.5];
nsamples = 100;
data = draw_from_ellipsoid(B, mu, nsamples);
x = data(:,1);
y = data(:,2);

% calculate volume of bounding ellipsoid
const = pi^(D/2)/gamma(D/2 + 1);
VS = const*sqrt(det(B));

% partition points into subclusters using optimal_ellipsoids
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

% display number of subclusters
fprintf('TEST 1: K = %d ellipses, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot(x,y,'b+');
axis equal
for k = 1:K
  [u, v, a, b, vol] = ellipsoid_2d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  plot(u,v,'r')
end
print -depsc2 test_optimal_ellipsoids1.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 2: two separated ellipsoids, 2d
clear B mu Bs mus x y z data
D = 2;

% bounding matrices and centroids for two separated ellipsoids
B1 = [1 -0.5; -0.5 1];
mu1 = [0 0];
B2 = [0.5 .25; .25 1];
mu2 = [4 3];

% calculate volumes of bounding ellipsoids
const = pi^(D/2)/gamma(D/2 + 1);
vol1 = const*sqrt(det(B1));
vol2 = const*sqrt(det(B2));

% generate uniform samples in each ellipsoid
nsamples1 = 100;
data1 = draw_from_ellipsoid(B1, mu1, nsamples1);
x1 = data1(:,1);
y1 = data1(:,2);

nsamples2 = round(nsamples1*vol2/vol1); % scale so uniform
data2 = draw_from_ellipsoid(B2, mu2, nsamples2);
x2 = data2(:,1);
y2 = data2(:,2);

% combine simulated points
x=[x1; x2];
y=[y1; y2];
data = [x y];

% partition points into subclusters using optimal_ellipsoids
VS = vol1 + vol2;
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

% display number of subclusters
fprintf('TEST 2: K = %d ellipses, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot(x,y,'b+');
axis equal
for k = 1:K
  [u, v, a, b, vol] = ellipsoid_2d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  plot(u,v,'r')
end
print -depsc2 test_optimal_ellipsoids2.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 3: circle of connected ellipsoids, 2d
clear B mu Bs mus x y z data
D = 2;

x = [];
y = [];
VS = 0;
N = 12;
% generate uniform samples from the ellipsoids
for ii=1:N
  r = 4;
  phi = 2*pi*(ii-1)/N;
  mu = [r*cos(phi) r*sin(phi)];
  B = [0.25 0; 0 1];
  R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
  B = R' * B * R;
  nsamples = 100;
  data = draw_from_ellipsoid(B, mu, nsamples);
  x = [x; data(:,1)];
  y = [y; data(:,2)];

  % calculate volume of bounding ellipsoids
  const = pi^(D/2)/gamma(D/2 + 1);
  vol = const*sqrt(det(B));
  VS = VS+vol;
end

data = [x y];

% partition points into subclusters using optimal_ellipsoids
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

% display number of subclusters
fprintf('TEST 3: K = %d ellipses, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot(x,y,'b+');
axis equal
for k = 1:K
  [u, v, a, b, vol] = ellipsoid_2d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  plot(u,v,'r')
end
print -depsc2 test_optimal_ellipsoids3.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 4: single ellipsoid, 3d
clear B mu Bs mus x y z data
D = 3;

% generate uniform samples from a single ellipsoid
B = [1 0.5 0; 0.5 1 -0.2; 0 -0.2 1];
mu = [1 2 3];
nsamples = 500;
data = draw_from_ellipsoid(B, mu, nsamples);
x = data(:,1);
y = data(:,2);
z = data(:,3);

% calculate volume of bounding ellipsoid
const = pi^(D/2)/gamma(D/2 + 1);
VS = const*sqrt(det(B));

% partition points into subclusters using optimal_ellipsoids
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

% display number of subclusters
fprintf('TEST 4: K = %d ellipsoids, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot3(x, y, z, 'b+');
alpha(.1);
axis equal
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
for k = 1:K
  [u, v, w, a, b, c, vol] = ellipsoid_3d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  colormap([1 0 0]);
  surf(u,v,w)
  alpha(.1) % almost completely transparent
  axis equal
end
print -depsc2 test_optimal_ellipsoids4.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 5: two separated ellipsoids, 3d
clear B mu Bs mus x y z data
D = 3;

% bounding matrices and centroids for two separated ellipsoids
B1 = [1 0.5 0; 0.5 1 -0.2; 0 -0.2 1];
mu1 = [1 2 3];
B2 = [0.5 0 0.2;  0 0.5 0; 0.2 0 .25];
mu2 = [0 0 0];

% calculate volumes of bounding ellipsoids
const = pi^(D/2)/gamma(D/2 + 1);
vol1 = const*sqrt(det(B1));
vol2 = const*sqrt(det(B2));

% generate uniform samples in each ellipsoid
nsamples1 = 500;
data1 = draw_from_ellipsoid(B1, mu1, nsamples1);
x1 = data1(:,1);
y1 = data1(:,2);
z1 = data1(:,3);

nsamples2 = round(nsamples1*vol2/vol1); % scale so uniform
data2 = draw_from_ellipsoid(B2, mu2, nsamples2);
x2 = data2(:,1);
y2 = data2(:,2);
z2 = data2(:,3);

% combine simulated points
x=[x1; x2];
y=[y1; y2];
z=[z1; z2];
data = [x y z];

% partition points into subclusters using optimal_ellipsoids
VS = vol1 + vol2;
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

% display number of subclusters
fprintf('TEST 5: K = %d ellipsoids, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot3(x, y, z, 'b+');
alpha(.1);
axis equal
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
for k = 1:K
  [u, v, w, a, b, c, vol] = ellipsoid_3d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  colormap([1 0 0]);
  surf(u,v,w)
  alpha(.1)
  axis equal
end
print -depsc2 test_optimal_ellipsoids5.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 6: torus of connected ellipsoids, 3d
clear B mu Bs mus x y z data
D = 3;

x = [];
y = [];
z = [];
VS = 0;
N = 12;
% generate uniform samples from the ellipsoids
for ii=1:N
  r = 4;
  phi = 2*pi*(ii-1)/N;
  mu = [r*cos(phi) r*sin(phi) 0];
  B = [0.25 0 0; 0 1 0; 0 0 0.25];
  R = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
  B = R' * B * R;
  nsamples = 500;
  data = draw_from_ellipsoid(B, mu, nsamples);
  x = [x; data(:,1)];
  y = [y; data(:,2)];
  z = [z; data(:,3)];

  % calculate volume of bounding ellipsoids
  const = pi^(D/2)/gamma(D/2 + 1);
  vol = const*sqrt(det(B));
  VS = VS+vol;
end

data = [x y z];

% partition points into subclusters using optimal_ellipsoids
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

fprintf('TEST 6: K = %d ellipsoids, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot3(x, y, z, 'b+');
alpha(.1);
axis equal
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
for k = 1:K
  [u, v, w, a, b, c, vol] = ellipsoid_3d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  colormap([1 0 0]);
  surf(u,v,w)
  alpha(.1)
  axis equal
end
print -depsc2 test_optimal_ellipsoids6.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST 7: uniformly distributed points in 3d torus
clear B mu Bs mus x y z data

% generate uniform samples in a 3d torus in xy plane
a = 2; % inner radius
b = 3; % outer radius
npts = 6000;
data = draw_from_torus(2, 3, npts);
x = data(:,1);
y = data(:,2);
z = data(:,3);

% calculate volume of torus
d = 0.5*(a+b);
R = 0.5*(b-a);
VS = 2*pi*d*pi*R^2;

% partition points into subclusters using optimal_ellipsoids
[Bs, mus, VEs, ns] = optimal_ellipsoids(data, VS);
K = length(ns);

fprintf('TEST 7: K = %d ellipsoids, F(S) = sum(VE)/VS = %f\n', K, sum(VEs)/VS)

% plot points and bounding ellipsoids for each subcluster
figure
plot3(x, y, z, 'b+');
alpha(.1);
axis equal
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')
for k = 1:K
  [u, v, w, a, b, c, vol] = ellipsoid_3d(Bs((k-1)*D+1:k*D,:), mus(k,:));
  hold on
  colormap([1 0 0]);
  surf(u,v,w)
  alpha(.1)
  axis equal
end
print -depsc2 test_optimal_ellipsoids7.eps

return


% script to test draw from ellipsoid
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% TEST 1: 2-d ellipsoid
% define test parameters
N = 100;
B = [1 .9; .9 1];
mu = [2 3];

% calculate bounding ellisoid for reference
[x, y, a, b, vol] = ellipsoid_2d(B, mu);

% get points
pnts = draw_from_ellipsoid(B, mu, N);

% plot bounding ellipse and points
figure
plot(x, y, 'r');
hold on
plot(pnts(:,1), pnts(:,2), '+');
axis equal
xlabel('X-axis')
ylabel('Y-axis')

print -depsc2 test_draw_from_ellipsoid1.eps

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% TEST 2: 3-d ellipsoid
% define test parameters
N = 1000;
B = [1 0.5 0; 0.5 1 -0.2; 0 -0.2 1];
mu = [1 2 3];

% calculate bounding ellisoid for reference
[x, y, z, a, b, c, vol] = ellipsoid_3d(B, mu);

% get points
pnts = draw_from_ellipsoid(B, mu, N);

% plot bounding ellipsoid and points
figure
colormap([1 0 0]); % red
surf(x, y, z)
alpha(.1) % almost completely transparent
axis equal
hold on
plot3(pnts(:,1), pnts(:,2), pnts(:,3), '+');
axis equal
xlabel('X-axis')
ylabel('Y-axis')
zlabel('Z-axis')

print -depsc2 test_draw_from_ellipsoid2.eps


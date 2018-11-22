function [x, y, a, b, vol] = ellipsoid_2d(B, mu)
%
% calculate (x,y) coordinates of the bounding 2-d ellipsoid defined
% by the bounding matrix B and centroid mu:
%
% Bounding ellipsoid = {u in R^2 | (u-mu)T B^{-1} (u-mu) <= 1}
%
% a, b: semi-major, semi-minor axes
% vol = pi*a*b (volume)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the eigenvectors/values of the bounding matrix B
[V, E] = eig(B);

% calculate the sqrt of the (diagonal) eigenvalue matrix
D = sqrt(E);

% extract semi-major and semi-minor axes
a = max(diag(D));
b = min(diag(D));
vol = pi*a*b;

% discrete values of phi
phi = linspace(0, 2*pi, 100);

% unit vector pointing in 'phi' direction
nx = cos(phi);
ny = sin(phi);

% transform unit vector to ellipsoid
for ii=1:length(phi)

     temp = V*D*[nx(ii) ny(ii)]';
     x(ii)=temp(1);
     y(ii)=temp(2);

end

% displace by mu
x = x + mu(1);
y = y + mu(2);

%figure
%plot(x,y,'r');
%axis equal
%xlabel('X-axis')
%ylabel('Y-axis')

return


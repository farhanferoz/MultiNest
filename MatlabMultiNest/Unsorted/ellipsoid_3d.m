function [x, y, z, a, b, c, vol] = ellipsoid_3d(B, mu)
%
% calculate (x,y,z) coordinates of the bounding 3-d ellipsoid defined
% by the bounding matrix B and centroid mu:
%
% Bounding ellipsoid = {u in R^3 | (u-mu)T B^{-1} (u-mu) <= 1}
%
% a, b, c: lengths of principle axes
% vol = (4/3)*pi*a*b*c (volume)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the eigenvectors/values of the bounding matrix B
[V, E] = eig(B);

% calculate sqrt of diagonal eigenvalue matrix
D = sqrt(E);

% extract lengths of principle axes
temp = sort(diag(D), 'descend');
a = temp(1);
b = temp(2);
c = temp(3);
vol = (4/3)*pi*a*b*c;

% discrete values of theta, phi
theta = linspace(0, pi, 20);
phi = linspace(0, 2*pi, 20);
thetaT = transpose(theta);

% unit vector points pointing in 'phi' and 'theta'
nx = sin(thetaT) * cos(phi);
ny = sin(thetaT) * sin(phi);
nz = cos(thetaT) * ones(size(phi));

% transform unit 2-sphere vectors to ellipsoid
for ii=1:length(theta)
  for jj=1:length(phi)

     temp = V*D*[nx(ii,jj) ny(ii,jj) nz(ii,jj)]';
     x(ii,jj)=temp(1);
     y(ii,jj)=temp(2);
     z(ii,jj)=temp(3);

  end
end

% displace by mu
x = x + mu(1);
y = y + mu(2);
z = z + mu(3);

%figure
%colormap([1 0 0]); % red
%surf(x,y,z)
%alpha(.1) % almost completely transparent
%axis equal
%xlabel('X-axis')
%ylabel('Y-axis')
%zlabel('Z-axis')

return

function pnts = draw_from_torus(a, b, N )

% function pnts = draw_from_torus(a, b,  N )
%
% This function draws points uniformly from an 3d torus,
% which lies in the xy-plane with center at (0,0,0) and
% inner radius a, outer radius b.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make sure that b>a
if b<=a
    error('DRAW FROM TORUS: b must be > a\n');
end

% radius of circular cross-section of torus
R = 0.5*(b-a);

% displacement of circular cross-section from origin
D = 0.5*(b+a);

% random unit vectors in 2-d
pt = randn(N, 2);
fac = sum(pt.^2, 2);
u = zeros(N,2);
for i=1:N
    u(i,:) = pt(i,:)/sqrt(fac(i));
end

% random scale factors for points in circular cross-section
% (unit radius)
rs = rand(N, 1);
fac = rs.^(1/2);

% generate points in torus
pnts = zeros(N,3);
for i=1:N
    % multiply unit vectors by random scale factors
    temp = fac(i)*u(i,:);

    % scale by radius of circular cross-section and
    % displace from origin
    temp = R*temp + [D 0];

    % rotate around z-direction by phi (uniform in [0, 2 pi])
    rho = temp(1);
    z = temp(2);
    phi = 2*pi*rand;

    pnts(i,:) = [rho*cos(phi) rho*sin(phi) z];

end

return

function plot_2d_livepoints_with_ellipses(livepoints, Bs, mus)
%
% plot 2-d livepoints in unit square and bounding ellipses
%
% NOTE:  mainly used for debug purposes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% extract number of ellipsoids and dimension
K = size(mus, 1);
ndims = size(mus, 2);

% plot live points and bounding ellipse for ndims=2 only
if ndims==2
  for k = 1:K
      [u, v, a, b, vol] = ellipsoid_2d(squeeze(Bs(k,:,:)), squeeze(mus(k,:)));
      plot(u,v,'r')
      hold on
  end

  plot(livepoints(:,1), livepoints(:,2), '+')
  xlim([-0.5 1.5])
  ylim([-0.5 1.5])
  drawnow
end

return


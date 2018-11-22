function histmat  = hist2(x, y, xedges, yedges)
% function histmat  = hist2(x, y, xedges, yedges)
%
% Extract 2D histogram data containing the number of events
% of [x , y] pairs that fall in each bin of the grid defined by
% xedges and yedges. The edges are vectors with monotonically
% non-decreasing values.
%
%EXAMPLE
%
% events = 1000000;
% x1 = sqrt(0.05)*randn(events,1)-0.5; x2 = sqrt(0.05)*randn(events,1)+0.5;
% y1 = sqrt(0.05)*randn(events,1)+0.5; y2 = sqrt(0.05)*randn(events,1)-0.5;
% x= [x1;x2]; y = [y1;y2];
%
%For linearly spaced edges:
% xedges = linspace(-1,1,64); yedges = linspace(-1,1,64);
% histmat = hist2(x, y, xedges, yedges);
% figure; pcolor(xedges,yedges,histmat'); colorbar ; axis square tight ;
%
%For nonlinearly spaced edges:
% xedges_ = logspace(0,log10(3),64)-2; yedges_ = linspace(-1,1,64);
% histmat_ = hist2(x, y, xedges_, yedges_);
% figure; pcolor(xedges_,yedges_,histmat_'); colorbar ; axis square tight ;

% University of Debrecen, PET Center/Laszlo Balkay 2006
% email: balkay@pet.dote.hu

% Copyright (c) 2009, Laszlo Balkay
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%    * Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in
%      the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if nargin ~= 4
    error ('The four input arguments are required!');
    return;
end
if any(size(x) ~= size(y))
    error ('The size of the two first input vectors should be same!');
    return;
end

[xn, xbin] = histc(x,xedges);
[yn, ybin] = histc(y,yedges);

%xbin, ybin zero for out of range values
% (see the help of histc) force this event to the
% first bins
xbin(find(xbin == 0)) = inf;
ybin(find(ybin == 0)) = inf;

xnbin = length(xedges);
ynbin = length(yedges);

if xnbin >= ynbin
    xy = ybin*(xnbin) + xbin;
      indexshift =  xnbin;
else
    xy = xbin*(ynbin) + ybin;
      indexshift =  ynbin;
end

%[xyuni, m, n] = unique(xy);
xyuni = unique(xy);
xyuni(end) = [];
hstres = histc(xy,xyuni);
clear xy;

histmat = zeros(ynbin,xnbin);
histmat(xyuni-indexshift) = hstres;


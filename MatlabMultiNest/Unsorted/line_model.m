function y = line_model(x, parnames, parvals)

% y = line_model(x, parnames, parvals)
%
% This function will return a line given by the equation y = mx + b, where
% m is the line's gradient and b is its y-intersept. The
% input parameters are:
%   x - the x values at which y will be calculated
%   parnames - a cell array containing the parameters names. These can be
%       in any order, but must include the following parameters:
%           {'m', 'b'}
%   parvals - a cell array containing the values of the parameters given in
%       parnames. These must be in the same order as in parnames.
%
%--------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
%--------------------------------------------------------------------------

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length!');
end

nparams = lpn;

% extract parameter values
for ii=1:nparams
  switch parnames{ii}
    case 'm'
      m = parvals{ii};
    case 'b'
      b = parvals{ii};
  end
end

% calculate y-values
y = m*x + b;

return

function y = sinusoid_model(t, parnames, parvals)

% function y = sinusoid_model(t, parnames, parvals)
%
% This function creates a sinusiod given by the equation
%     y = amp * sin(2*pi*f*t + phi)
% The input parameters are:
%
% t - the discrete times (in sec) at which y will be calculated
% parnames - a cell array containing the parameter names.
%     These can be in any order, but must include the following
%         {'amp', 'phi', 'f'}
%     where amp is the amplitude, phi is the initial phase
%     (in radians), and f is the frequency (in Hz)
% parvals - a cell array containing the values of the parameters
%     given in parnames.  These must be in the same order as in
%     parnames
%
%--------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
%--------------------------------------------------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    case 'amp'
      amp = parvals{ii};
    case 'phi'
      phi = parvals{ii};
    case 'f'
      f = parvals{ii};
  end
end

% calculate sinusiod
y = amp * sin(2*pi*f*t + phi);

end

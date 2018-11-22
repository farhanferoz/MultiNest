function logL = logL_mixture_gaussian(data, model, parnames, parvals)

% logL = logL_mixture_gaussian(data, model, parnames, parvals)
%
% This function will compute the log likelihood of a mixture gaussian:
%
%   L = prod_{i=1}^N ( (1-Pb)/sqrt{(2*pi) sigma_i^2}
%           exp[ -(y_i - model(x_i,params))^2 / 2*sigma_i^2 ]
%         + Pb/sqrt{(2*pi)(Vb + sigma_i^2)}
%           exp[ -(y_i - Yb)^2 / 2*(Vb + sigma_i^2) ] )
%
% The input parameters are:
%   data - a cell array with three columns
%       1) x values, 2) y values, 3) sigma_y values
%   model - the function handle for the signal model
%   parnames - a cell array listing the names of the model parameters
%              and parameters (Pb, Yb, logVb) for the mixture gaussian
%   parvals - a cell array containing the values of the parameters given
%       in parnames. These must be in the same order as in parnames.
%       If parvals is an empty vector the noise-only likelihood will be
%       calculated.
%
% -------------------------------------------------------------------------
%           This is the format required by nested_sampler.m.
% -------------------------------------------------------------------------

% check whether model is a string or function handle
if ischar(model)
    fmodel = str2func(model);
elseif isa(model, 'function_handle')
    fmodel = model;
else
    error('Error... Expecting a model function!');
end

% get data values from cell array
x = cell2mat(data(1));
y = cell2mat(data(2));
s = cell2mat(data(3));
N = length(x);

% evaluate the model
if isempty(parvals)
    % if parvals is not defined get the null likelihood (noise model
    % likelihood)
    md = 0;
else
    md = feval(fmodel, x, parnames, parvals);

    % if the model returns a NaN then set the likelihood to be zero (e.g.
    % loglikelihood to be -inf
    if isnan(md)
        logL = -inf;
        return;
    end
end

% extract parameters for mixture gaussian
nparams = length(parnames);
for ii=1:nparams
  switch parnames{ii}
    case 'Pb'
      Pb = parvals{ii};
    case 'Yb'
      Yb = parvals{ii};
    case 'logVb'
      logVb = parvals{ii};
      Vb = exp(logVb);
  end
end

% calculate the log likelihood
logL = 0;
for ii=1:N
  logL = logL + ...
         log( (1-Pb)/sqrt(2*pi*s(ii)^2) * ...
              exp(-0.5*(y(ii)-md(ii))^2 / s(ii)^2) + ...
              Pb/sqrt(2*pi*(Vb+s(ii)^2)) * ...
              exp(-0.5*(y(ii)-Yb)^2 / (Vb+s(ii)^2)) );
end

if isnan(logL)
    error('Error: log likelihood is NaN!');
end

return

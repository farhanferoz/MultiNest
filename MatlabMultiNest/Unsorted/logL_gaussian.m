function logL = logL_gaussian(data, model, parnames, parvals)

% logL = logL_gaussian(data, model, parnames, parvals)
%
% This function will compute the log likelihood of a multivariate
% gaussian:
%
%     L = 1/sqrt((2 pi)^N det C)
%         exp[-0.5*(y - model(x,params))^T * inv(C) * (y - model(x,params))]
%
% The input parameters are:
%     data - a cell array with three columns
%            { x values, y values, C: covariance matrix }
%     NOTE: if C is a single number, convert to a diag covariance matrix
%     model - the function handle for the signal model.
%     parnames - a cell array listing the names of the model parameters
%     parvals - a cell array containing the values of the parameters given
%         in parnames. These must be in the same order as in parnames.
%         If parvals is an empty vector the noise-only likelihood will be
%         calculated.
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
x = data{1};
y = data{2};
C = data{3};
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

% get inverse of covariance matrix and the log of the determinant
if size(C,1)==1 || N == 1
    % variance is just a single value
    invC = 1/C;
    lDetC = N*log(C);
    % calculate the log likelihood
    logL = -0.5*(y - md)' * invC * (y - md);
elseif size(C,1) == N && size(C,2) == N && N > 1
    % use trick from http://xcorr.net/2008/06/11/log-determinant-of-positive-definite-matrices-in-matlab/
    % to calculate log of determinant and avoid infinities
    Cchol = chol(C); % Cholesky decomposition
    lDetC = 2*sum(log(diag(Cchol)));

    % calculate the log likelihood (use / for matrix inverse)
    logL = -0.5*((y - md)'/C)*(y - md);
end

% calculate the log likelihood
logL = logL - 0.5*N*log(2*pi) - 0.5*lDetC;

if isnan(logL)
    error('Error: log likelihood is NaN!');
end

return

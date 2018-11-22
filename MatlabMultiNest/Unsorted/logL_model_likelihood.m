function logL = logL_model_likelihood(data, model, parnames, parvals)

% check whether model is a string or function handle
if ischar(model)
    fmodel = str2func(model);
elseif isa(model, 'function_handle')
    fmodel = model;
else
    error('Error... Expecting a model function!');
end

x = data;

md = feval(fmodel, x, parnames, parvals);

logL = log(md);
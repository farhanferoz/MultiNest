function y = triangle_model(~, parnames, parvals)

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
    case 'gradient'
      grad = parvals{ii};
    case 'x'
      x = parvals{ii};
  end
end

y = x*grad;
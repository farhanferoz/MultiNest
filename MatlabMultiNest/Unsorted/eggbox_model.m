function y = eggbox_model(~, parnames, parvals)

% An "eggbox"-like model with N dimensions

% check that parnames and parvals have the same length
lpn = length(parnames);
lpv = length(parvals);
if lpn ~= lpv
    error('Error: parnames and parvals are not the same length!');
end

nparams = lpn;

pars = zeros(nparams,1);

% extract parameter values
for ii=1:nparams
  pars(ii) = parvals{ii};
end

chi = 1;
for i=1:nparams
	x = pars(i)*10.0*pi;
    chi = chi*cos(x/2.0);
end

y = exp((chi + 2.0)^5.0);
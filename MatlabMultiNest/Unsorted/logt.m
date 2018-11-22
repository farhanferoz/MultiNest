function y = logt(N)
%
% This function returns the logarithm y=log(t) of the largest of
% N independent samples drawn from a unif distribution in (0,1)
%
% P(t) = N t^{N-1} has E(log t) = -1/N, dev(log t) = 1/N
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = max(rand(N,1));
y = log(t);

return


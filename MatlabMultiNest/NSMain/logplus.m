function logz  = logplus(logx, logy)
%
% logz = logplus(logx, logy)
%
% Given logx and logy, this function returns logz=log(x+y).
% It avoids problems of dynamic range when the
% exponentiated values of x or y are very large or small.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isinf(logx) && isinf(logy)
    logz = -inf;
    return;
end

if logx > logy
    logz = logx+log(1.+exp(logy-logx));
else
    logz = logy+log(1.+exp(logx-logy));
end

return


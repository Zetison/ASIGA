function x = normalize_lp(x)
if nargin < 2
    p = 2;
end
x = x./vecnorm(x,p,1);
x(isnan(x)) = 0;
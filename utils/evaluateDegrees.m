function P = evaluateDegrees(xi,p)

P = zeros(p+1,1);

for i = 1:p+1
    P(i) = xi^(i-1);
end
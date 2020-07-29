function I = findArcLength(R_a,f,a,b)

fun = @(theta) sqrt(R_a^2 - f^2*cos(theta).^2);

I = integral(fun,a,b);
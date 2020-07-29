
%-------------------------------------------------------------------------------
function [x_hat,iter]=golden_ratio(f,x_min,x_max,tol)
r1 = (sqrt(5)-1)/2;
r2 = r1^2;

x = zeros(1,4);

x(1) = x_min;
x(4) = x_max;
x(2) = x(1) + (x(4)-x(1))*r2;
x(3) = x(1) + (x(4)-x(1))*r1;

F = arrayfun(f,x);

iter = 0;
while x(4)-x(1) > tol,
    iter = iter+1;
    if F(2) > F(3),
        x=[x(1),x(1)+(x(3)-x(1))*r2,x(2),x(3)];
        F=[F(1),f(x(2)),F(2),F(3)];
    else
        x=[x(2),x(3),x(2)+(x(4)-x(2))*r1,x(4)];
        F=[F(2),F(3),f(x(3)),F(4)];
    end
end
x_hat = x(2);

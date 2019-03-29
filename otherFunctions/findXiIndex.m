function xi = findXiIndex(x_fixed,nurbs,epsilon)

xi_min = 0;
xi_max = 1;

xi = 0.5;
v = evaluateNURBS(nurbs, xi);
x = v(1);

while abs(x_fixed - x) > epsilon
    if x > x_fixed
        xi_max = xi;
    else
        xi_min = xi;
    end
    xi = (xi_max+xi_min)/2;
    v = evaluateNURBS(nurbs, xi);
    x = v(1);
end
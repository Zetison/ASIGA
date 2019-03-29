function r = bisection(f, a, b, N, eps_step)
% Assume that f(a) and f(b) has different sign, and tries to find a root in between
if f(a) == 0
    f_a = f(a+eps_step);
else
    f_a = f(a);
end
if f(b) == 0
    f_b = f(b-eps_step);
else
    f_b = f(b);
end
if f_a*f_b > 0
    error('f(a) and f(b) should have different sign')
end
% We will iterate N times and if a root was not
% found after N iterations, an exception will be thrown.
for k = 1:N
    % Find the mid-point
    c = (a + b)/2;
    f_c = f(c);
    if f_c == 0
        r = c;
        return;
    elseif f_c*f_a < 0
        b = c;
    else
        a = c;
    end

    if b - a < eps_step
        if abs(f_a) < abs(f_b)
            r = a;
            return;
        else
            r = b;
            return;
        end
    end
end

error('The method did not converge');
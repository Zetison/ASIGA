function [x,i] = newtonsMethod(f,dfdx,x_0,n,Eps,bnd)
% This algorithm assumes x = O(1)
x_prev = x_0;
f_x_prev = f(x_0);
if nargin < 6
    for i = 1:n
        dfdx_x = dfdx(x_prev);
        if abs(dfdx_x) < Eps
            dfdx_x = dfdx(x_prev+Eps);
        end
        f_x = f(x_prev);
        x = x_prev - f_x/dfdx_x;
        if abs(f_x) > abs(f_x_prev)
            break
        end
        if abs(x-x_prev) < Eps
            return
        end
        x_prev = x;
        f_x_prev = f_x;
    end
else
    for i = 1:n
        dfdx_x = dfdx(x_prev);
        if abs(dfdx_x) < Eps
            x_prev = bnd(1)+rand()*bnd(2);
            dfdx_x = dfdx(x_prev);
        end
        f_x = f(x_prev);
        x = x_prev - f_x/dfdx_x;
        if x < bnd(1)
            x = bnd(1);
        end
        if x > bnd(2)
            x = bnd(2);
        end
        if abs(x-x_prev) < Eps
            return
        end
        x_prev = x;
    end
end
warning('Newtons method did not converge')
% error('Newtons method did not converge')
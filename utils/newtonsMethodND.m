function [xi,i] = newtonsMethodND(f,dfdx,xi_0,n,Eps,bnd)
xi = xi_0;
d_p = numel(xi_0);
f_xi = evaluateFunctions(f,dfdx,xi);
d = numel(f_xi);
for i = 1:n
    [f_xi,dfdx_xi] = evaluateFunctions(f,dfdx,xi);
    while any(isnan(dfdx_xi(:))) || any(isinf(dfdx_xi(:))) || cond(dfdx_xi) > 1/Eps
        xi = bnd(:,1)+rand(d_p,1).*(bnd(:,2)-bnd(:,1));
        [f_xi,dfdx_xi] = evaluateFunctions(f,dfdx,xi);
    end
    xi_prev = xi;
    f_xi_prev = f_xi;
    if d_p == 1 && d > 1 % use Moore-Penrose pseudoinverse
        xi = xi - dfdx_xi.'/(dfdx_xi.'*dfdx_xi)*f_xi;
    else
        xi = xi - dfdx_xi\f_xi;
    end
    for j = 1:d_p
        if xi(j) < bnd(j,1)
            xi(j) = bnd(j,1);
        elseif xi(j) > bnd(j,2)
            xi(j) = bnd(j,2);
        end
    end
            
    if norm(xi-xi_prev)/norm(xi) < Eps
        return
    end    
    if norm(f_xi) > norm(f_xi_prev)
        warning('A function value got largerNewtons method did not converge')
        break
    end
end
warning('Newtons method did not converge')

% error('Newtons method did not converge')

function [f_xi,dfdx_xi] = evaluateFunctions(f,dfdx,xi)
if isnan(dfdx)
    F = f(xi);
    f_xi = F{1};
    dfdx_xi = F{2};
else
    f_xi = f(xi);
    dfdx_xi = dfdx(xi);
end
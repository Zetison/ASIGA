function [x,i] = newtonsMethodND(f,dfdx,x_0,n,Eps,bnd)
x_prev = x_0;
% dfdx_prev = Inf;
for i = 1:n
    dfdx_x = dfdx(x_prev);
    if any(any(isnan(dfdx_x))) || any(any(isinf(dfdx_x))) || cond(dfdx_x) > 1/Eps
%         warning('Newtons method did not converge')
        x_prev(1) = bnd(1,1)+rand()*(bnd(1,2)-bnd(1,1));
        x_prev(2) = bnd(2,1)+rand()*(bnd(2,2)-bnd(2,1));
        dfdx_x = dfdx(x_prev);
%         x = x_0;
%         return
    end
    x = x_prev - dfdx_x\f(x_prev);
    for j = 1:numel(x_0)
        if x(j) < bnd(j,1)
            x(j) = bnd(j,1);
        elseif x(j) > bnd(j,2)
            x(j) = bnd(j,2);
        end
    end
            
    if norm(x-x_prev)/norm(x) < Eps
        return
    end    
%     if norm(dfdx_x) > norm(dfdx_prev)
%         return
%     end
    x_prev = x;
%     dfdx_prev = dfdx_x;
end
% warning('Newtons method did not converge')
% error('Newtons method did not converge')
function [xi,eta] = newtonsMethod2D(F,J,x_0,n,Eps,bnd)
xi = x_0(:,1);
eta = x_0(:,2);
xi_prev = xi;
eta_prev = eta;
compact = ~iscell(J);
% dfdx_prev = Inf;
for i = 1:n
%     dfdx_x = dfdx(x_prev);
%     if any(any(isnan(dfdx_x))) || any(any(isinf(dfdx_x))) || cond(dfdx_x) > 1/Eps
% %         warning('Newtons method did not converge')
%         x_prev(1) = bnd(1,1)+rand()*(bnd(1,2)-bnd(1,1));
%         x_prev(2) = bnd(2,1)+rand()*(bnd(2,2)-bnd(2,1));
%         dfdx_x = dfdx(x_prev);
% %         x = x_0;
% %         return
%     end
    if compact
        FF = F(xi,eta);
        F1 = FF{1};
        F2 = FF{2};
        J11 = FF{3};
        J12 = FF{4};
        J22 = FF{5};
    else
        J11 = J{1}(xi,eta);
        J12 = J{2}(xi,eta);
        J22 = J{3}(xi,eta);
        F1 = F{1}(xi,eta);
        F2 = F{2}(xi,eta);
    end
    D = J11.*J22 - J12.^2;
    xi = xi_prev - (J22.*F1 - J12.*F2)./D;
    eta = eta_prev - (-J12.*F1 + J11.*F2)./D;
    if nargin == 6
        xi(xi < bnd(1,1)) = bnd(1,1);
        xi(xi > bnd(1,2)) = bnd(1,2);
        eta(eta < bnd(2,1)) = bnd(2,1);
        eta(eta > bnd(2,2)) = bnd(2,2);
    end
            
    if (norm(xi-xi_prev)/norm(xi)+norm(eta-eta_prev)/norm(eta))/2 < Eps
        return
    end    
%     if norm(dfdx_x) > norm(dfdx_prev)
%         return
%     end
    xi_prev = xi;
    eta_prev = eta;
%     dfdx_prev = dfdx_x;
end
% warning('Newtons method did not converge')
% error('Newtons method did not converge')
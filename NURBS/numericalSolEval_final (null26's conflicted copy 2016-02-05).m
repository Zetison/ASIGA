function [u, v, dudx, dudy, dudz, J] = numericalSolEval_final(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights, controlPts, U)
maxJacobian = 1e4;
if nargout <= 2
    R = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
else
    [R, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
end


n = length(Xi) - (p+1);
m = length(Eta) - (q+1);
l = length(Zeta) - (r+1);

i1 = findKnotSpan(n, p, xi, Xi);
i2 = findKnotSpan(m, q, eta, Eta);
i3 = findKnotSpan(l, r, zeta, Zeta);

u = zeros(1,size(U,2));
v = zeros(1,3);
pts = zeros(length(R),3);
counter = 1;
for k3 = 1:r+1
    A3 = i3 - r + k3 - 1;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            A = (m*n)*(A3-1) + n*(A2-1) + A1;   
            u = u + U(A,:)*R(counter);
            v = v + controlPts(A,:)*R(counter);
            pts(counter,:) = controlPts(A,:);
            counter = counter + 1;
        end
    end
end
u = u.';
v = v';
% if eta == 0.5 && zeta == 0
%     keyboard
% end
if nargout > 2
    J  = pts' * [dRdxi' dRdeta' dRdzeta']; 
    epsilon = 1e-15;
    
    xi_orig = xi;
    eta_orig = eta;
    zeta_orig = zeta;
    while cond(J) > maxJacobian
%         if xi_orig == 0
%             xi = xi_orig + epsilon;
%         else
%             xi = xi_orig - epsilon;
%         end
%         if eta_orig == 0
%             eta = eta_orig + epsilon;
%         else
%             eta = eta_orig - epsilon;
%         end
        if zeta_orig == 0
            zeta = zeta_orig + epsilon;
        else
            zeta = zeta_orig - epsilon;
        end
        [~, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
        epsilon = epsilon*10;
    end           
    dRdX = J'\[dRdxi; dRdeta; dRdzeta];
    
    dudx = zeros(1,size(U,2));
    dudy = zeros(1,size(U,2));
    dudz = zeros(1,size(U,2));
    
    counter = 1;
    for k3 = 1:r+1
        A3 = i3 - r + k3 - 1;
        for k2 = 1:q+1
            A2 = i2 - q + k2 - 1;
            for k1 = 1:p+1
                A1 = i1 - p + k1 - 1;
                A = (m*n)*(A3-1) + n*(A2-1) + A1;   

                dudx = dudx + dRdX(1,counter)*U(A,:);
                dudy = dudy + dRdX(2,counter)*U(A,:);
                dudz = dudz + dRdX(3,counter)*U(A,:);
                counter = counter + 1;
            end
        end
    end
end





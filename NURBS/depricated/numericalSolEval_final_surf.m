function [u, v, dudx, dudy, J] = numericalSolEval_final_surf(xi, eta, p, q, Xi, Eta, weights, controlPts, U)
error('Depricated. Use evalNURBSsol() instead')
maxJacobian = 1e6;
if nargout <= 2
    R = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
else
    [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
end

u = zeros(1,size(U,2));
v = zeros(1,3);
pts = zeros(length(R),3);
counter = 1;
for k2 = 1:q+1
    for k1 = 1:p+1
        u = u + U(counter,:)*R(counter);
        v = v + controlPts(counter,:)*R(counter);
        pts(counter,:) = controlPts(counter,:);
        counter = counter + 1;
    end
end
v = v';

if nargout > 2
    J  = pts' * [dRdxi' dRdeta']; 
    epsilon = 1e-15;
    
    xi_orig = xi;
    eta_orig = eta;
    while cond(J) > maxJacobian
        if xi_orig == 0
            xi = xi_orig + epsilon;
        else
            xi = xi_orig - epsilon;
        end
        if eta_orig == 0
            eta = eta_orig + epsilon;
        else
            eta = eta_orig - epsilon;
        end
        if xi < 0 || xi > 1 || eta < 0 || eta > 1
            warning('A singularity of the geometric mapping was not properly handled in Post-Processing.')
            break
        end
        [~, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
        J  = pts' * [dRdxi' dRdeta'];
        epsilon = epsilon*10;
    end           
    dRdX = J'\[dRdxi; dRdeta];
    dudx = zeros(1,size(U,2));
    dudy = zeros(1,size(U,2));
    
    counter = 1;
    for k2 = 1:q+1
        for k1 = 1:p+1
            dudx = dudx + dRdX(1,counter)*U(counter,:);
            dudy = dudy + dRdX(2,counter)*U(counter,:);
            counter = counter + 1;
        end
    end
end





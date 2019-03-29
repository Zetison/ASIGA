function [u, v, dudx, dudy] = numericalSolEval_2D_final(xi, eta, p, q, Xi, Eta, weights, controlPts, U, computeDeriv)
maxJacobian = 1e5;
if nargin < 10
    computeDeriv = false;
    R = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
else
    [R, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
end


n = length(Xi) - (p+1);
m = length(Eta) - (q+1);

i1 = findKnotSpan(n, p, xi, Xi);
i2 = findKnotSpan(m, q, eta, Eta);

u = zeros(1,size(U,2));
v = zeros(1,2);
pts = zeros(length(R),2);
counter = 1;
for k2 = 1:q+1
    A2 = i2 - q + k2 - 1;
    for k1 = 1:p+1
        A1 = i1 - p + k1 - 1;
        A = n*(A2-1) + A1;   
        u = u + U(A,:)*R(counter);
        v = v + controlPts(A,:)*R(counter);
        pts(counter,:) = controlPts(A,:);
        counter = counter + 1;
    end
end
u = u.';
v = v';

if computeDeriv
    J  = pts' * [dRdxi' dRdeta']; 
    if cond(J) > maxJacobian
        if xi == 0
            xi = xi + 1e-2;
            [~, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
            J  = pts' * [dRdxi' dRdeta'];
        else
            xi = xi - 1e-2;
            [~, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
            J  = pts' * [dRdxi' dRdeta'];
        end
    end
    if cond(J) > maxJacobian
        if eta == 0
            eta = eta + 1e-2;
            [~, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
            J  = pts' * [dRdxi' dRdeta'];
        else
            eta = eta - 1e-2;
            [~, dRdxi, dRdeta] = NURBS2DBasis(xi, eta, p, q, Xi, Eta, weights);
            J  = pts' * [dRdxi' dRdeta'];
        end
    end      
    dRdX = J'\[dRdxi; dRdeta];
    
    dudx = zeros(1,2);
    dudy = zeros(1,2);
    
    counter = 1;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            A = (m*n)*(A3-1) + n*(A2-1) + A1;   

            dudx = dudx + dRdX(1,counter)*U(A,:);
            dudy = dudy + dRdX(2,counter)*U(A,:);
            counter = counter + 1;
        end
    end
end





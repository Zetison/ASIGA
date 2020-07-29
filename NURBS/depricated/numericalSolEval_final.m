function [u, v, dudx, dudy, dudz, J] = numericalSolEval_final(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights, controlPts, U)
error('Depricated. Use evalNURBSsol() instead')

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
            A = m*n*(A3-1) + n*(A2-1) + A1;   
            u = u + U(A,:)*R(counter);
            v = v + controlPts(A,:)*R(counter);
            pts(counter,:) = controlPts(A,:);
            counter = counter + 1;
        end
    end
end
v = v';

if nargout > 2
%     symbolicPrecision = (abs(eta) < 10*eps || abs(eta-1) < 10*eps);
    dudx = zeros(1,size(U,2));
    dudy = zeros(1,size(U,2));
    dudz = zeros(1,size(U,2));
        
    J  = pts' * [dRdxi' dRdeta' dRdzeta'];
    symbolicPrecision = cond(J) > 1e10;
%     symbolicPrecision =  false;
    if cond(J) > 1e10 && ~symbolicPrecision
        return
    end
    if symbolicPrecision
        digits(100)
        xi = vpa(xi); eta = vpa(eta); zeta = vpa(zeta); 
        Xi = vpa(Xi); Eta = vpa(Eta); Zeta = vpa(Zeta);
        weights = vpa(weights);
        pts = vpa(pts);
        
        if abs(xi) < 10*eps
            xi = xi + 1e5*eps;
        elseif abs(xi-1) < 10*eps
            xi = xi - 1e5*eps;
        end
        if abs(eta) < 10*eps
            eta = eta + 1e5*eps;
        elseif abs(eta-1) < 10*eps
            eta = eta - 1e5*eps;
        end
        if abs(zeta) < 10*eps
            zeta = zeta + 1e5*eps;
        elseif abs(zeta-1) < 10*eps
            zeta = zeta - 1e5*eps;
        end
        [~, dRdxi, dRdeta, dRdzeta] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);
        J  = pts' * [dRdxi' dRdeta' dRdzeta'];
    end
    dRdX = J'\[dRdxi; dRdeta; dRdzeta];
    
    if symbolicPrecision
        dudx = vpa(dudx);
        dudy = vpa(dudy);
        dudz = vpa(dudz);
        U = vpa(U);
    end
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





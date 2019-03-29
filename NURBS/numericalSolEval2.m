function [u, v] = numericalSolEval2(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights, Ux, Uy, Uz, P)

[R, ~, ~, ~] = NURBS3DBasis(xi, eta, zeta, p, q, r, Xi, Eta, Zeta, weights);

n = length(Xi) - (p+1);
m = length(Eta) - (q+1);
l = length(Zeta) - (r+1);

i1 = findKnotSpan(n, p, xi, Xi);
i2 = findKnotSpan(m, q, eta, Eta);
i3 = findKnotSpan(l, r, zeta, Zeta);

u = zeros(3,1);
counter = 1;
for k3 = 1:r+1
    A3 = i3 - r + k3 - 1;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            A = (m*n)*(A3-1) + n*(A2-1) + A1;   
            u = u + [Ux(A); Uy(A); Uz(A)]*R(counter);
            counter = counter + 1;
        end
    end
end

v = zeros(3,1);
counter = 1;
for k3 = 1:r+1
    A3 = i3 - r + k3 - 1;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            v = v + R(counter)*P(1:3, A1, A2, A3);
            counter = counter + 1;
        end
    end
end
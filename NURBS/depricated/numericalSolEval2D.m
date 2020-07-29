function [u, v, dudx, dudy, d2udx2, d2udy2] = numericalSolEval2D(xi, eta, varCol, U)

error('Depricated. Use evalNURBSsol() instead')
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
weights = varCol.weights;
controlPts = varCol.controlPts;

[R, dRdxi, dRdeta, d2Rdxi2, d2Rdeta2] = NURBS2DBasis2(xi, eta, p, q, Xi, Eta, weights);

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
        u = u + U(A)*R(counter);
        v = v + controlPts(A,:)*R(counter);
        pts(counter,:) = controlPts(A,:);
        counter = counter + 1;
    end
end
v = v';

if nargout > 2
    J  = pts' * [dRdxi' dRdeta']; 
    
    dRdX = J'\[dRdxi; dRdeta];
    %% This is not general
    d2RdX2 = [d2Rdxi2/J(1,1)^2; d2Rdeta2/J(2,2)^2];
    dudx = 0;
    dudy = 0;
    d2udx2 = 0;
    d2udy2 = 0;
    
    counter = 1;
    for k2 = 1:q+1
        A2 = i2 - q + k2 - 1;
        for k1 = 1:p+1
            A1 = i1 - p + k1 - 1;
            A = n*(A2-1) + A1;   

            dudx = dudx + dRdX(1,counter)*U(A);
            dudy = dudy + dRdX(2,counter)*U(A);
            d2udx2 = d2udx2 + d2RdX2(1,counter)*U(A);
            d2udy2 = d2udy2 + d2RdX2(2,counter)*U(A);
            counter = counter + 1;
        end
    end
end





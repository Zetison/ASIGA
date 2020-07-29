function [u, v, dudx, d2udx2, d3udx3] = numericalSolEval1D(xi, varCol, U)

error('Depricated. Use evalNURBSsol() instead')
p = varCol.nurbs.degree;
Xi = varCol.nurbs.knots;
weights = varCol.weights;
controlPts = varCol.controlPts;
[R, dRdxi, d2Rdxi2, d3Rdxi3] = NURBS1DBasis2(xi, p, Xi, weights);


n = length(Xi) - (p+1);

i1 = findKnotSpan(n, p, xi, Xi);
pts = zeros(length(R),1);

u = 0;
v = 0;
counter = 1;
for k1 = 1:p+1
    A = i1 - p + k1 - 1;
    u = u + U(A)*R(counter);
    v = v + controlPts(A)*R(counter);
    pts(counter) = controlPts(A,:);
    counter = counter + 1;
end

if nargout > 2
    J  = dRdxi*pts; 
    J2  = d2Rdxi2*pts; 
    J3  = d3Rdxi3*pts; 
    dRdX = dRdxi/J;
    d2RdX2 = d2Rdxi2/J^2 - dRdxi*J2/J^3;
    d3RdX3 = d3Rdxi3/J^3 + 2*d2Rdxi2*J2/J - d2Rdxi2*J2/J^4 - dRdxi*(J3*J - 3*J2^2)/J^5;
    dudx = 0;
    d2udx2 = 0;
    d3udx3 = 0;
    
    counter = 1;
    for k1 = 1:p+1
        A = i1 - p + k1 - 1;

        dudx = dudx + dRdX(counter)*U(A);
        d2udx2 = d2udx2 + d2RdX2(counter)*U(A);
        d3udx3 = d3udx3 + d3RdX3(counter)*U(A);
        counter = counter + 1;
    end
end





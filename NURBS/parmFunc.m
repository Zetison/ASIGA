function nurbs = parmFunc(Xi,p,f)

xi_arr = aveknt(Xi, p+1);

n = length(Xi)-(p+1);
A = zeros(n);
dummy = f(0);
F = zeros(n,numel(dummy));
for i = 1:n
    xi = xi_arr(i);
    i1 = findKnotSpan(n, p, xi, Xi);
    A(i,i1-p:i1) = BsplineBasis(i1, xi, p, Xi, 0);
    F(i,:) = f(xi);
end
coeffs = (A\F).';
coeffs(end+1,:) = 1;
nurbs = createNURBSobject(coeffs, Xi);
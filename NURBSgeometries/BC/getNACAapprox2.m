function cntrlPts = getNACAapprox2(t,p,Xi)

f = @(xi) getNACA(xi,t);

xi_arr = aveknt(Xi, p+1);

n = length(Xi)-(p+1);
A = zeros(n);
F = zeros(n,2);
for i = 1:n
    xi = xi_arr(i);
    i1 = findKnotSpan(n, p, xi, Xi);
    A(i,i1-p:i1) = Bspline_basis(i1, xi, p, Xi, 0);
    F(i,:) = [xi^2, f(xi)];
end

[xi_rep, I] = getRepeatedKnots(Xi,p);
% if idx
%     
% else
    cntrlPts = A\F;
% end
cntrlPts(I,:) = [(xi_rep.^2).', f(xi_rep).'];
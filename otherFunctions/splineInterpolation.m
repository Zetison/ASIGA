function splineStrct = splineInterpolation(f,xmin,xmax,n,p)



Xi = [zeros(1,p) linspace(0,1,n-p+1) ones(1,p)];

idx = 2:n-1;
if p == 1
    xi_arr = linspace(0,1,n);
elseif p == 2
    xi_arr = [0 (Xi(idx+2).*(Xi(idx+3)-Xi(idx+1))+Xi(idx+1).*(Xi(idx+2)-Xi(idx)))./(Xi(idx+3)-Xi(idx+1)+Xi(idx+2)-Xi(idx)) 1];
else
    xi_arr = zeros(1,n);
    xi_arr(end) = 1;
    for i = 2:n-1
        ximin = Xi(i);
        ximax = Xi(i+p+1);
        B = @(xi) evalBspline(i,xi,n,p,Xi);
        xi_arr(i) = golden_ratio(B,ximin+eps,ximax-eps,eps);
    end
end
    
F = arrayfun(f,xmin+xi_arr*(xmax-xmin))';
 
A = sparse(n,n);

for i = 1:n
    xi = xi_arr(i);
    i1 = findKnotSpan(n, p, xi, Xi);
    N = Bspline_basis(i1, xi, p, Xi, 0);
    A(i,i1-p:i1) = N;
end
coeffs = A\F;
coeffs = [coeffs'; ones(1,n)];
splineStrct = createNURBSobject(coeffs,Xi);

function B = evalBspline(i,xi,n,p,Xi)

i1 = findKnotSpan(n, p, xi, Xi);
B = Bspline_basis(i1, xi, p, Xi, 0);
B = B(p+1+i-i1);


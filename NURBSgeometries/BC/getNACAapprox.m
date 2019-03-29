function [cntrlPts, Xi] = getNACAapprox(t,p,N,xi_T,f_T)

f = @(xi) getNACA(xi,t);

if nargin < 3
    N = [1,p+1];
    xi_T = [0,1];
    f_T = [0,0];
end
Xi = zeros(1,p+1);
for i = 2:numel(N)
    xi = xi_T(i);
    Xi = [Xi, linspace2(xi_T(i-1),xi,N(i)-N(i-1)-p), xi*ones(1,p)];
end
Xi = [Xi,1];

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

if numel(N) > 3
    cntrlPts = zeros(n,2);
    cntrlPts(:,1) = A\F(:,1);
    i = 3;
    y = f_T(i);
    idx = [N(i)-1,N(i)+1];
    F(:,2) = F(:,2) - A(:,idx)*ones(2,1)*y;
    F(idx,:) = [];
    A(idx,:) = [];
    A(:,idx) = [];
    temp = A\F(:,2);
    cntrlPts(setdiff(1:n,idx),2) = temp;
    cntrlPts(idx,2) = ones(2,1)*y;
else
    cntrlPts = A\F;
end
cntrlPts(2,1) = 0; % for exactness
cntrlPts(N,:) = [xi_T.^2.',f_T.']; % for exactness
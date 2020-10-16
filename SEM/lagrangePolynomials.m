function [p,dpdx] = lagrangePolynomials(xi,jj,N,Xi)
if isrow(xi)
    xi = xi';
end
p = zeros(numel(xi),numel(jj));
for i = 1:numel(jj)
    j = jj(i);
    indices = [1:j-1,j+1:N];
    p(:,i) = prod((xi-Xi(indices))./repmat(Xi(j)-Xi(indices),numel(xi),1),2);
end
if nargout > 1
    dpdx = zeros(numel(xi),numel(jj));
    for i = 1:numel(jj)
        j = jj(i);
        dpdx(:,i) = lagrangePolynomialsNthDeriv(xi,j,N,Xi,1);
    end
end
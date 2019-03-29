function p = lagrangePolynomials(xi,j,N,Xi)
if isrow(xi)
    xi = xi';
end
indices = [1:j-1,j+1:N];
p = prod((xi-Xi(indices))./repmat(Xi(j)-Xi(indices),numel(xi),1),2);
function dpdxi = lagrangePolynomialsDeriv(xi,i,N,Xi)
if isrow(xi)
    xi = xi';
end

dpdxi = zeros(size(xi));
for j = 1:N
    if j ~= i
        if i < j
            l = [1:i-1,i+1:j-1,j+1:N];
        else
            l = [1:j-1,j+1:i-1,i+1:N];
        end
        dpdxi = dpdxi + prod((xi-Xi(l))./repmat(Xi(i)-Xi(l),numel(xi),1),2)/(Xi(i)-Xi(j));
    end
end
function nrm = norm2(x)
% x is a MxN matrix. Computes the norms nrm(i) = norm(x(i,:))

nrm = sqrt(sum(abs(x).^2,2));
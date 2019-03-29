function nrm = norm3(x)
% x is a MxN matrix. Computes the norms nrm(i) = norm(x(i,:))

M = size(x,1);

nrm = zeros(M,1,class(x));
for i = 1:M
    nrm(i) = sqrt(x(i,:)*x(i,:).');
end
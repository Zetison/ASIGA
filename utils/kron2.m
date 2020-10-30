function C = kron2(A,B)
% Assumes size(A,1) == size(B,1)
[m,n1] = size(A);
n2 = size(B,2);
C = zeros(m*n2,n1);
for i = 1:n1
    C(:,i) = reshape((A(:,i).*B).',m*n2,1);
end
C = reshape(C.',n1*n2,m);
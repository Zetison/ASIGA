function y = matrixArrayProd(A,x)
% y = matrixArrayProd(A,x) computes the product y(:,i,j,...) = A*x(:,i,j,...)
sizes = size(x);

y = reshape(A*x(:,:),[size(A,1),sizes(2:end)]);
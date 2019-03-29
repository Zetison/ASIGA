function P = elementProd(v, e)
% v is a Nx1 matrix and e is a NxM matrx. Computes the products P(:,i) =
% v.*e(:,i). If e is a 1xM matrix, the kronecker product kron(v,e) is
% computed
if size(e,1) == 1
    P = kron(v,e);
else
    P = zeros(size(e),class(v));
    M = size(e,2);
    for i = 1:M
        P(:,i) = v.*e(:,i);
    end
end


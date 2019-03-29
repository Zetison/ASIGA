function B = dot3(A,n)
% If A is a MxN matrix and n is a Nx1 vector. Computes all products
% B(i) = A(i,:)*n. That is, dot products without complex conjugation.
% Similar case if A is a 1xN vector and n is a NxM matrix


M = size(A,1);
if M > 1
    if size(n,2) > 1
        B = zeros(M,1,class(A));
        for i = 1:M
            B(i) = A(i,:)*n(:,i);
        end
    else
        B = zeros(M,1,class(A));

        for i = 1:M
            B(i) = A(i,:)*n;
        end
    end
else
    M = size(n,2);
    B = zeros(M,1,class(A));

    for i = 1:M
        B(i) = A*n(:,i);
    end
end
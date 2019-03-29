function x = blockMatrixSolver(A,b,blockIndices)

switch length(blockIndices)
    case 0
        x = A\b;
    case 1
        m = blockIndices(1);
        A11_in = inv(A(1:m,1:m));
        A21_11 = A(m+1:end,1:m)/A(1:m,1:m);
        x2 = (A(m+1:end,m+1:end) - A21_11*A(1:m,m+1:end))\(b(m+1:end)-A21_11*b(1:m));
        x1 = A(1:m,1:m)\(b(1:m)-A(1:m,m+1:end)*x2);
        x = [x1; x2];
    otherwise
        error('Not implemented')
end
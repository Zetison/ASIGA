function B = copyVector(A,n,type)

m = length(A);
B = zeros(m*n,1);

switch type
    case 1
        for i = 1:n
            idx = (1:m) + m*(i-1);
            B(idx) = A;
        end
    case 2
        for i = 1:m
            idx = (1:n) + n*(i-1);
            B(idx) = A(i)*ones(n,1);
        end
end

        
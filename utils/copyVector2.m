function B = copyVector2(A,n1,n2,type)

m = length(A);
B = zeros(m*n1*n2,1);

switch type
    case 1
        B = repmat(A.',n1*n2,1);
    case 2
        B = reshape(repmat(A,n1,1,n2),m*n1*n2,1);
    case 3
        B = reshape(repmat(A,n1*n2,1),m*n1*n2,1);
end

        
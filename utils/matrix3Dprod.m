function p = matrix3Dprod(A,B)

p = zeros(size(A,1),size(B,2),size(B,3));

for i = 1:size(B,2)
    for j = 1:size(B,3)
        p(:,i,j) = A*B(:,i,j);
    end
end
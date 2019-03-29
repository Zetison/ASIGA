function B = uniquetol2(A,tol)

B = zeros(size(A));
counter = 1;
for i = 1:size(A,1)
    foundElement = false;
    for j = 1:i-1
        if abs(A(i,:) - A(j,:)) < tol
            foundElement = true;
            break
        end
    end
    if ~foundElement
        B(counter,:) = A(i,:);
        counter = counter + 1;
    end
end
B = B(1:counter-1,:);
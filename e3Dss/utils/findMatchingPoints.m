function [indices_x, indices_y] = findMatchingPoints(x,y,Eps)


M = size(x,1);
N = size(y,1);
indices_x = zeros(max(M,N),1);
indices_y = zeros(max(M,N),1);
counter = 1;
X = kron(x,ones(N,1));
Y = kron(ones(M,1),y);
normXY = norm2(X-Y) < 10*Eps;
for i = 1:M
    for j = 1:N
        if normXY(j+(i-1)*N)
            indices_x(counter) = i;
            indices_y(counter) = j;
            counter = counter + 1;
        end
    end
end
indices_x = indices_x(1:counter-1);
indices_y = indices_y(1:counter-1);
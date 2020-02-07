function xi = aveknt(Xi,k)

n = numel(Xi)-k;
xi = zeros(1,n);
p = k-1;
j = 1:n;
for i = 1:p
    xi = xi+Xi(j+i);
end
xi = xi/p;
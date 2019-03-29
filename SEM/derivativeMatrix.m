function D = derivativeMatrix(xi,N)

D = zeros(N);
for i = 1:N
    for alpha = 1:N
        if alpha ~= i
            if i < alpha
                l = [1:i-1,i+1:alpha-1,alpha+1:N];
            else
                l = [1:alpha-1,alpha+1:i-1,i+1:N];
            end
            D(i,alpha) = prod((xi(alpha)-xi(l))./(xi(i)-xi(l)))/(xi(i)-xi(alpha));
        end
    end
end
D(1,1) = -N*(N-1)/4;
D(N,N) = N*(N-1)/4;

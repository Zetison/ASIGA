function B = adaptiveChebychevInterp(f,epsilon)
N = 2;
B = zeros(N+1,1);
F = zeros(N+1,1);
for j = 0:N
    x_j = cos(j*pi/N);
    F(j+1) = f(x_j);
end
    
for k = 0:N
    temp = 0;
    for j = 0:N
        x_j = cos(j*pi/N);
        if j == 0 || j == N
            temp = temp + 0.5*chebyshevTilde(k,x_j)*F(j+1);
        else
            temp = temp + chebyshevTilde(k,x_j)*F(j+1);
        end
    end
    B(k+1) = 2/N*temp;
end

while abs(B(end)) > epsilon*sum(abs(B))
    N = 2*N;
    Ftemp = F;
    B = zeros(N+1,1);
    F = zeros(N+1,1);
    F(1:2:end) = Ftemp;
    for j = 1:2:N-1
        x_j = cos(j*pi/N);
        F(j+1) = f(x_j);
    end
    for k = 0:N
        temp = 0;
        for j = 0:N
            x_j = cos(j*pi/N);
            if j == 0 || j == N
                temp = temp + 0.5*chebyshevTilde(k,x_j)*F(j+1);
            else
                temp = temp + chebyshevTilde(k,x_j)*F(j+1);
            end
        end
        B(k+1) = 2/N*temp;
    end
end
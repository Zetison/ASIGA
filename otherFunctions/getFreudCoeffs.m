function [B, S] = getFreudCoeffs(n, r)

a = @(m) gamma((m+1)/(r+1))/(r+1);

B = zeros(n+1,1);
if isa(r,'sym')
    B = vpa(B);
end
B(1) = 1;

A = zeros(n,n+1);
if isa(r,'sym')
    A = vpa(A);
end

for i = 0:n-1
    for j = 0:n
        A(i+1,j+1) = (-1)^j*a(i+j);
    end
end
if n == 0
    B(1) = 1;
else
    c = det(A(:,2:end));
    for i = 1:n
        B(i+1) = det(A(:,[1:i,i+2:end]))/c;
    end
end

if nargout == 2
    S = 0;
    for i = 0:2*n
        temp = 0;
        if i <= n
            for j = 0:i
                temp = temp + B(j+1)*B(i-j+1);
            end
        else
            for j = 0:2*n-i
                temp = temp + B(n-j+1)*B(i+j-n+1);
            end
        end
        S = S + a(i)*temp;
    end    
end

function z = findAllZeros(f,range,n)

x = linspace(range(1), range(2), n);
z = [];
counter = 0;
for i = 1:n-1
    if sign(f(x(i))) ~= sign(f(x(i+1)))
        z(counter+1) = fzero(f,(x(i)+x(i+1))/2);
        counter = counter + 1;
    end
end
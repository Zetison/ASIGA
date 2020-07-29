function A = nonUniformLinspace(a,b,n,B)

if isempty(B)
    A = linspace(a,b,n);
    return
end
if a == B(1,1)
    A = a;
else
    A = linspace(a,B(1,1),n);
end

for i = 1:size(B,1)
    if i ~= size(B,1)
        A = [A linspace2(B(i,1), B(i,2), B(i,3)) linspace(B(i,2), B(i+1,1),n)];  
    else
        if B(end,2) == b
            A = [A linspace2(B(i,1), B(i,2), B(i,3)) b];  
        else
            A = [A linspace2(B(i,1), B(i,2), B(i,3)) linspace(B(i,2), b,n)];  
        end
    end
end



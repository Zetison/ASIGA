function u = addAdjacentValues(v,epsilon)

u = zeros(2*length(v),1);
for i = 1:length(v)
    u(2*i-1) = v(i) - epsilon;
    u(2*i) = v(i) + epsilon;
end
    
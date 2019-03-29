function nrm = norm4(x)

sz = size(x);

if length(sz) == 2
    nrm = zeros(size(x,1),1);
    for i = 1:size(x,1)
        nrm(i) = norm(x(i,:));
    end
end
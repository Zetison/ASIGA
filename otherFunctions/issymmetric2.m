function [issymmetr, idx] = issymmetric2(A,newEpsilon)
idx = zeros(1,2);
if issparse(A)
    [idxRows,idxColumns,values] = find(A);
    maxA = max(abs(values));
    issymmetr = true;
    for i = 1:length(values)
        if abs(values(i) - A(idxColumns(i), idxRows(i)))/maxA > newEpsilon
            issymmetr = false;
            idx = [idxRows(i), idxColumns(i)];
            return
        end
    end
        
else
    maxA = max(max(abs(A)));
    N = size(A,1);
    issymmetr = true;
    for i = 1:N
        for j = i:N
            if abs(A(i,j) - A(j,i))/maxA > newEpsilon
                issymmetr = false;
                idx = [i, j];
                return
            end
        end
    end
end
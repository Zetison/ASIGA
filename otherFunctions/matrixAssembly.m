function A = matrixAssembly(values, idx, n_en, noDofs, noElems, type)
switch type
    case 1
        A = zeros(n_en, noDofs);
        for e = 1:noElems
            for j = 1:n_en
                A(:,idx(j,e)) = A(:,idx(j,e)) + values(:,j,e);
            end
        end
    case 2
        A = zeros(noDofs);
        for e = 1:noElems
            for i = 1:n_en
                A(idx(i,e),:) = A(idx(i,e),:) + values(i,:,e);
            end
        end        
end

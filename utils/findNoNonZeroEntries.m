function noNonZeroEntries = findNoNonZeroEntries(noDofs,element,noElems,dimension)

K_temp = repmat(struct('rs',cell(1)),noDofs,1);
M_temp = repmat(struct('rs',cell(1)),noDofs,1);
for e = 1:noElems
    sctr = element(e,:);          %  element scatter vector
    n_en = length(sctr);
    for i = 1:n_en
        ii = sctr(i);
        for j = 1:n_en
            jj = sctr(j);
            if isempty(K_temp(ii).rs)
                K_temp(ii).rs(1) = jj;
                M_temp(ii).rs(1) = jj;
            else
                elementExist = false;
                for k = 1:length(K_temp(ii).rs)
                    if K_temp(ii).rs(k) == jj
                        elementExist = true;
                        break;
                    end
                end
                if ~elementExist
                    K_temp(ii).rs(k+1) = jj;
                    M_temp(ii).rs(k+1) = jj;
                end
            end
        end
    end
end
noNonZeroEntries = 0;
for i = 1:noDofs
    noNonZeroEntries = noNonZeroEntries + length(K_temp(i).rs);
end

noNonZeroEntries = dimension*noNonZeroEntries;
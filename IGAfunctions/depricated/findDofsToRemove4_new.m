function varCol = findDofsToRemove4_new(varCol, homDirichletDofs)
error('Depricated, use findDofsToRemove instead')

if length(varCol.nurbs.knots) == 3
    n = varCol.nurbs.number(1);
    m = varCol.nurbs.number(2);
    l = varCol.nurbs.number(3);
    childrenNodes = varCol.childrenNodes;

    counter = 1;
    counter2 = 1;
    gluedNodes = cell(1,2*l + (m-2)*l);
    for k = 1:l
        for j = 1:m
            for i = 1:n
                if j == 1 || j == m
                    gluedNodes{counter2} = counter:(counter+n-1);
                    counter2 = counter2 + 1;
                    counter = counter + n;
                    break
                elseif i == 1
                    gluedNodes{counter2} = [counter (counter+n-1)];   
                    counter2 = counter2 + 1;
                end
                counter = counter + 1;
            end
        end
    end
else
    n = varCol.nurbs.number(1);
    m = varCol.nurbs.number(2);
    childrenNodes = varCol.childrenNodes;

    counter = 1;
    counter2 = 1;
    gluedNodes = cell(1,m);
    for j = 1:m
        for i = 1:n
            if j == 1 || j == m
                gluedNodes{counter2} = counter:(counter+n-1);
                counter2 = counter2 + 1;
                counter = counter + n;
                break
            elseif i == 1
                gluedNodes{counter2} = [counter (counter+n-1)];   
                counter2 = counter2 + 1;
            end
            counter = counter + 1;
        end
    end
end

d = varCol.dimension;
dofsToRemove = zeros(1,length(childrenNodes)*d);
for i = 1:d
    dofsToRemove(i:d:end) = d*(childrenNodes-1)+i;
end
if ~isempty(homDirichletDofs)
    dofsToRemove = [dofsToRemove homDirichletDofs];
end

dofsToRemove = sort(unique(dofsToRemove));

varCol.dofsToRemove = dofsToRemove;
varCol.gluedNodes = gluedNodes;
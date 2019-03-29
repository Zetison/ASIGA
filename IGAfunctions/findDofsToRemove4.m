function [dofsToRemove, gluedNodes] = findDofsToRemove4(noCtrlPts, homDirichletDofs, n,m,l, dimension,childrenNodes)

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


dofsToRemove = [];
for d = 1:dimension
    dofsToRemove = [dofsToRemove childrenNodes+(d-1)*noCtrlPts];
end
if ~isempty(homDirichletDofs)
    dofsToRemove = [dofsToRemove homDirichletDofs];
end

dofsToRemove = sort(unique(dofsToRemove));
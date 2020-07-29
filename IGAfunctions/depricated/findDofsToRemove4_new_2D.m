function varCol = findDofsToRemove4_new_2D(varCol, homDirichletDofs)
error('Depricated, use findDofsToRemove instead')
noCtrlPts = varCol.noCtrlPts;
n = varCol.nurbs.number(1);
m = varCol.nurbs.number(2);
childrenNodes = varCol.childrenNodes;

counter = 1;
counter2 = 1;
gluedNodes = cell(1,m);
for j = 1:m
    for i = 1:n
        if i == 1
            gluedNodes{counter2} = [counter (counter+n-1)];   
            counter2 = counter2 + 1;
        end
        counter = counter + 1;
    end
end


dofsToRemove = [];
for d = 1:varCol.dimension
    dofsToRemove = [dofsToRemove childrenNodes+(d-1)*noCtrlPts];
end
if ~isempty(homDirichletDofs)
    dofsToRemove = [dofsToRemove homDirichletDofs];
end

dofsToRemove = sort(unique(dofsToRemove));

varCol.dofsToRemove = dofsToRemove;
varCol.gluedNodes = gluedNodes;
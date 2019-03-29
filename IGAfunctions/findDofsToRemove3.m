function [element, dofsToRemove, gluedNodes] = findDofsToRemove3(element, controlPts, noCtrlPts, newEpsilon, homDirichletDofs, dimension,maxConnections,maxConnected, candidatesIndices)

gluedNodes = findMatchIndices3(controlPts, newEpsilon,maxConnections,maxConnected,candidatesIndices);

childrenNodes = zeros(1,length(cell2mat(gluedNodes)) - length(gluedNodes));

counter = 1;
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        indices = (element == gluedNodes{i}(j));
        element(indices) = parentIdx;
        childrenNodes(counter) = gluedNodes{i}(j);
        counter = counter + 1;
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
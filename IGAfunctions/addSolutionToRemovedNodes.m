function U = addSolutionToRemovedNodes(U_temp,dofsToRemove,noDofs,gluedNodes,noCtrPts,homDirichletDofs)

nonDirichletNodes = setdiff(1:noDofs,dofsToRemove');

U = zeros(noDofs,1);

U(nonDirichletNodes) = U_temp;

for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        U(gluedNodes{i}(j)) = U(parentIdx);
        U(gluedNodes{i}(j)+noCtrPts) = U(parentIdx+noCtrPts);
        U(gluedNodes{i}(j)+2*noCtrPts) = U(parentIdx+2*noCtrPts);
    end
end

if ~isempty(homDirichletDofs)
    U(homDirichletDofs) = 0;
end
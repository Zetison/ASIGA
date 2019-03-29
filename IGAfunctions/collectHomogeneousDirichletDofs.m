function homDirichletDofs = collectHomogeneousDirichletDofs(rigidNodes, xConstNodes, yConstNodes, zConstNodes, noCtrlPts)


xConstNodes = [xConstNodes rigidNodes];
yConstNodes = [yConstNodes rigidNodes];
zConstNodes = [zConstNodes rigidNodes];

xConstNodes = unique(xConstNodes);
yConstNodes = unique(yConstNodes);
zConstNodes = unique(zConstNodes);

homDirichletDofs = [];
if ~isempty(xConstNodes)
    homDirichletDofs = [homDirichletDofs xConstNodes];
end
if ~isempty(yConstNodes)
    homDirichletDofs = [homDirichletDofs yConstNodes+noCtrlPts];
end
if ~isempty(zConstNodes)
    homDirichletDofs = [homDirichletDofs zConstNodes+2*noCtrlPts];
end
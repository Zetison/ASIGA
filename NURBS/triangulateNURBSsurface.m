function [TRI,ppPts,patches,dofsToRemovePP] = triangulateNURBSsurface(patches,extraPts)

for patch_i = 1:numel(patches)
    patch = patches{patch_i};
    nurbs = patch.nurbs;
    
    tempXi = insertUniform(nurbs.knots{1},extraPts);
    tempEta = insertUniform(nurbs.knots{2},extraPts);
    nxi = numel(tempXi);
    neta = numel(tempEta);
    nPts = nxi*neta;
    prmPts = [tempXi, zeros(nxi,1); 
              ones(neta-1,1), tempEta(2:end);
              tempXi(end-1:-1:1), ones(nxi-1,1);
              zeros(neta-2,1), tempEta(end-1:-1:2);
              copyVector(tempXi(2:end-1),neta-2,1), copyVector(tempEta(2:end-1),nxi-2,2)];
          
    Eps = 1e5*eps;
    
    X = zeros(size(prmPts,1),3);
    for gp = 1:nPts
        X(gp,:) = evaluateNURBS(nurbs, prmPts(gp,:));
    end
    
    [~, I] = uniquetol(X,Eps,'ByRows',true, 'DataScale',max(norm2(X)));
    I = sort(I);
    X = X(I,:);
    prmPts = prmPts(I,:);
    noBndryPts = 2*nxi+2*neta-4-(nxi*neta-numel(I));
    
    nNormalPts = 5;
    normals = zeros(nNormalPts^2,3);
    [~,Q2D] = gaussianQuadNURBS(nNormalPts,nNormalPts); 
    Q2D = (1+Q2D)/2; 
    for gp = 1:size(Q2D,1)
        [~, dXdxi, dXdeta] = evaluateNURBS_deriv(nurbs, Q2D(gp,:));
        crossProd = cross(dXdxi, dXdeta); 
        normals(gp,:) = crossProd/norm(crossProd);
    end 
    n_avg = sum(normals,1)/nNormalPts^2;
    X_m = orthogonalTransform(X, n_avg);
    
    DT = delaunayTriangulation(X_m(:,1:2),[(1:noBndryPts).', ([2:noBndryPts,1]).']);
    IO = isInterior(DT);
    TRI = DT.ConnectivityList(IO,:);
    
    patches{patch_i}.ppPrmPts = prmPts;
    patches{patch_i}.ppPts = X;
    patches{patch_i}.noPts = size(X,1);
    patches{patch_i}.TRI = TRI;
    patches{patch_i}.noTRI = size(TRI,1);
%     if patch_i == 140
%         keyboard
%     end
    if any(norm2(X_m(:,1:2)-DT.Points) > 1e-12)
        error('delaunay has shifted the ordering of the points')
    end
%     plotNURBS(nurbs,{'alphaValue',0.8});
%     trisurf(TRI,X(:,1),X(:,2),X(:,3))
%     axis equal
%     hold on
%     ax = gca;               % get the current axis
%     ax.Clipping = 'off';    % turn clipping off
%     axis off
end


noPatches = numel(patches);
noTRIPatch = zeros(noPatches,1);
noPtsPatch = zeros(noPatches,1);
for i = 1:noPatches
    noTRIPatch(i) = patches{i}.noTRI;
    noPtsPatch(i) = patches{i}.noPts;
end


noCtrlPts = sum(noPtsPatch);
noElems = sum(noTRIPatch);
TRI = zeros(noElems,3);
ppPts = zeros(noCtrlPts,size(patches{1}.controlPts,2));
e = 1;
jC = 1;
maxDof = 0;
for i = 1:noPatches
    TRI(e:e+noTRIPatch(i)-1,:) = maxDof + patches{i}.TRI;
    ppPts(jC:jC+noPtsPatch(i)-1,:) = patches{i}.ppPts;
    maxDof = maxDof + noPtsPatch(i);
    e = e + noTRIPatch(i);
    jC = jC + noPtsPatch(i);
end

[~, gluedNodes] = uniquetol(ppPts,Eps,'ByRows',true, 'DataScale',max(norm2(ppPts)), 'OutputAllIndices', true);
repeatedNode = zeros(numel(gluedNodes),1);
for i = 1:numel(gluedNodes)
    repeatedNode(i) = numel(gluedNodes{i}) - 1;
end
gluedNodes(repeatedNode == 0) = [];
noChildrenNodes = sum(repeatedNode);

childrenNodes = zeros(1,noChildrenNodes);
counter = 1;
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        childrenIdx = gluedNodes{i}(j);
        indices = (TRI == childrenIdx);
        TRI(indices) = parentIdx;
        childrenNodes(counter) = childrenIdx;        
        counter = counter + 1;
    end
end

dofsToRemovePP = childrenNodes;
dofsToRemovePP = sort(unique(dofsToRemovePP));

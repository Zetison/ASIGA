function nurbs = cleanNURBS(nurbs,parms)

parms = [parms, 0, 1, 1/sqrt(2)];
cutOff = 1e6;
Eps = cutOff*eps;
for i = 1:numel(nurbs)
    for j = 1:numel(parms)
        nurbs{i}.coeffs(abs(nurbs{i}.coeffs - parms(j)) < Eps) = parms(j);
        nurbs{i}.coeffs(abs(nurbs{i}.coeffs + parms(j)) < Eps) = -parms(j);
    end
end

noPatches = numel(nurbs);
noCtrlPtsPatch = zeros(noPatches,1);
for i = 1:noPatches
    nPts = prod(nurbs{i}.number);
    noCtrlPtsPatch(i) = nPts;
end


noCtrlPts = sum(noCtrlPtsPatch);
controlPts = zeros(noCtrlPts,3);
jC = 1;
for i = 1:noPatches
    nPts = prod(nurbs{i}.number);
    controlPts(jC:jC+noCtrlPtsPatch(i)-1,:) = reshape(nurbs{i}.coeffs(1:3,:,:),3,nPts).';
    jC = jC + noCtrlPtsPatch(i);
end
[~, gluedNodes] = uniquetol(controlPts,Eps,'ByRows',true, 'DataScale',max(norm2(controlPts)), 'OutputAllIndices', true);
repeatedNode = zeros(numel(gluedNodes),1);
for i = 1:numel(gluedNodes)
    repeatedNode(i) = numel(gluedNodes{i}) - 1;
end
gluedNodes(repeatedNode == 0) = [];

for i = 1:length(gluedNodes)
    controlPts(gluedNodes{i}(2:end)) = controlPts(gluedNodes{i}(1));
end

counter = 1;
for i = 1:numel(nurbs)
    nPts = prod(nurbs{i}.number);
    nurbs{i}.coeffs(1:3,:,:) = reshape(controlPts(counter:counter+nPts-1,:).',3,nurbs{i}.number(1),nurbs{i}.number(2));
    counter = counter + nPts;
end
    
        
        
        
        
        
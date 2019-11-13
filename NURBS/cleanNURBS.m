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

d = size(nurbs{1}.coeffs,1)-1;

noCtrlPts = sum(noCtrlPtsPatch);
controlPts = zeros(noCtrlPts,d);
jC = 1;
for i = 1:noPatches
    nPts = prod(nurbs{i}.number);
    controlPts(jC:jC+noCtrlPtsPatch(i)-1,:) = reshape(nurbs{i}.coeffs(1:d,:,:),d,nPts).';
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
    nurbs{i}.coeffs(1:d,:,:,:) = reshape(controlPts(counter:counter+nPts-1,:).',[d,nurbs{i}.number]);
    counter = counter + nPts;
end
    
        
        
        
        
        
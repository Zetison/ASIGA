function nurbs = cleanNURBS(nurbs,parms,Eps)

parms = [parms, 0, 1, 0.5, 0.25, 0.125, 1/sqrt(2)];
if nargin < 3
    Eps = 1e6*eps;
end
for i = 1:numel(nurbs)
    for j = 1:numel(parms)
        nurbs{i}.coeffs(abs(nurbs{i}.coeffs - parms(j)) < Eps) = parms(j);
        nurbs{i}.coeffs(abs(nurbs{i}.coeffs + parms(j)) < Eps) = -parms(j);
        for ii = 1:numel(nurbs{i}.knots)
            nurbs{i}.knots{ii}(abs(nurbs{i}.knots{ii} - parms(j)) < Eps) = parms(j);
            nurbs{i}.knots{ii}(abs(nurbs{i}.knots{ii} + parms(j)) < Eps) = -parms(j);
        end
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
    controlPts(gluedNodes{i}(2:end),:) = repmat(controlPts(gluedNodes{i}(1),:),numel(gluedNodes{i})-1,1);
end

counter = 1;
for i = 1:numel(nurbs)
    nPts = prod(nurbs{i}.number);
    nurbs{i}.coeffs(1:d,:,:,:) = reshape(controlPts(counter:counter+nPts-1,:).',[d,nurbs{i}.number]);
    counter = counter + nPts;
end
    
        
        
        
        
        
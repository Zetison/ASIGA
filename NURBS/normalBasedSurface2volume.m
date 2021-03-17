function [nurbsVol,nurbs2] = normalBasedSurface2volume(nurbs,D)
if ~iscell(nurbs)
    nurbs = {nurbs};
end
task.misc.method = 'BA';
task.misc.progressBars = false;
task.misc.extraGP = [0,0,0];
task.varCol{1}.dimension = 3;
task.varCol{1}.nurbs = nurbs;
task.varCol{1} = findDofsToRemove(generateIGAmesh(convertNURBS(task.varCol{1})));
task.varCol{1}.operator = 'laplace';
task.varCol{1}.applyBodyLoading = true;
task.varCol{1}.buildMassMatrix = true;
task.varCol{1}.buildStiffnessMatrix = false;
task.varCol{1}.fieldDimension = 3;
task.varCol{1}.force = @(X,n) X + D*n;
task.varCol{1}.extraGP = 3*ones(1,2);
task.varCol{1}.progressBars = false;
task = buildMatrices(task,1);
task.varCol{1}.media = 'fluid';
[task,FF,~,~,A2] = collectMatrices(task);
UU = A2\FF;
task.varCol = postProcessSolution(task.varCol,UU);
noPatches = numel(nurbs);
nurbs2 = cell(1,noPatches);
counter = 1;
for i = 1:noPatches
    sz = size(nurbs{i}.coeffs);
    coeffs = zeros(sz);
    coeffs(4,:,:) = nurbs{i}.coeffs(4,:,:);
    noDofs = task.varCol{1}.noCtrlPtsPatch(i)*task.varCol{1}.dimension;
    coeffs(1:3,:,:) = reshape(task.varCol{1}.U(counter:counter+noDofs-1,:).',3,sz(2),sz(3));
    counter = counter + noDofs;
    nurbs2(i) = createNURBSobject(coeffs,nurbs{i}.knots);
end
nurbsVol = loftNURBS({nurbs,nurbs2});
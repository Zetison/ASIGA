% function nurbs = meshOptimizer %(nurbs,rigidNodes)
clear all
close all

study = 7;
m = 2;      % objective function exponent
objectiveFun = 3;
useJacobian = 1;
L = 2;
d_p = 3;
d = d_p;
nurbs = getPrismData('L',L,'d_p', d_p, 'd', d);
s = 0.75*L;
plotAffinity = 0;
plotOrthogonality = 1;
plotMeshQuality = 0;
Eps = 1e-10;
maxIter = 1e5;
extraGP = [3,3,3];

tol = 1e-10;
TolX = 1e-10;
TolFun = tol^2;
plotIterations = false;
prePlot = 1;
postPlot = true;
startMatlabPool

resolution = 32*ones(1,3);
% resolution = 100*ones(1,3);
% resolution = 2*ones(1,3);
% method = 'Nelder-Mead';
method = 'Gradient decent';
Lvec = zeros(1,d);
Lvec(1) = L;
if d_p == 3
    switch study
        case 1 % Simple linear box with only one CP to be optimized
            nurbs2 = nurbs;
            nurbs2{1}.coeffs(1:3,end,end,end) = [s,s,s];
        case 2 % Two boxes with irregularities at one patch to be optimized
            nurbs = elevateNURBSdegree(nurbs,1);
            nurbs2 = nurbs;
            nurbs2{1}.coeffs([1,3],end,end-1,end) = [s,s];
            nurbs2{1}.coeffs([2,3],end-1,end,end) = [s,s];
            nurbs2{1}.coeffs([1,2],end,end,end-1) = [s,s];
        case 3 % Same as case 2 but with a unit cube added
            nurbs = elevateNURBSdegree(nurbs,1);
            nurbs = [nurbs,nurbs];
            nurbs(2) = translateNURBS(nurbs(2),[L,0,0]);
            nurbs2 = nurbs;
            nurbs2{2}.coeffs([1,3],end,end-1,end) = [L+s,s];
            nurbs2{2}.coeffs([2,3],end-1,end,end) = [s,s];
            nurbs2{2}.coeffs([1,2],end,end,end-1) = [L+s,s];
        case 4 % Two boxes with irregularities at one patch to be optimized (including a shared CP)
            nurbs = elevateNURBSdegree(nurbs,1);
            nurbs = [nurbs,nurbs];
            nurbs(2) = translateNURBS(nurbs(1),[L,0,0]);
            nurbs2 = nurbs;
            nurbs2{2}.coeffs([1,3],end,end-1,end) = [L+s,s];
            nurbs2{2}.coeffs([2,3],end-1,end,end) = [s,s];
            nurbs2{2}.coeffs([1,2],end,end,end-1) = [L+s,s];
            nurbs2{2}.coeffs(1:3,1,end,end) = [1.75,1.75,1.75];
            nurbs2{1}.coeffs(1:3,end,end,end) = nurbs2{2}.coeffs(1:3,1,end,end);
        case 5 % Two boxes with irregularities at one patch to be optimized (including a shared CP)
            nurbs = [nurbs,nurbs];
            nurbs(2) = translateNURBS(nurbs(2),[L,0,0]);
            nurbs2 = nurbs;
            nurbs2{2}.coeffs(1:3,end,end,end) = [L+s,s,s];
            nurbs2{2}.coeffs(1:3,end-1,end,end) = [1.75,1.75,1.75];
            nurbs2{1}.coeffs(1:3,end,end,end) = nurbs2{2}.coeffs(1:3,end-1,end,end);
        case 6 % Load from file
            load('test.mat')
            nurbs2 = nurbs;
        case 7 % Load from file
            load('test_S1_vol.mat')
            nurbs2 = nurbs;
    end
else
    switch study
        case 1 % Simple linear box with only one CP to be optimized
            nurbs2 = nurbs;
            nurbs2{1}.coeffs(1:2,end,end) = [s,s];
        case 2 % Two boxes with irregularities at one patch to be optimized
            nurbs = elevateNURBSdegree(nurbs,1);
            nurbs2 = nurbs;
            nurbs2{1}.coeffs(1:2,end,end-1) = [1.1*s,s];
            nurbs2{1}.coeffs(1:2,end-1,end) = [s,1.1*s];
        case 3 % Same as case 2 but with a unit square added
            nurbs = elevateNURBSdegree(nurbs,1);
            nurbs = [nurbs,nurbs];
            nurbs(2) = translateNURBS(nurbs(2),Lvec);
            nurbs2 = nurbs;
            nurbs2{2}.coeffs(1:2,end,end-1) = [L+1.1*s,s];
            nurbs2{2}.coeffs(1:2,end-1,end) = [L+s,1.1*s];
        case 4 % Two boxes with irregularities at one patch to be optimized (including a shared CP)
            nurbs = elevateNURBSdegree(nurbs,1);
            nurbs = [nurbs,nurbs];  
            nurbs(2) = translateNURBS(nurbs(1),Lvec);
            nurbs2 = nurbs;
            nurbs2{2}.coeffs(1:2,end,end-1) = [L+1.1*s,s];
            nurbs2{2}.coeffs(1:2,end-1,end) = [L+s,1.1*s];
            nurbs2{2}.coeffs(1:2,1,end) = [1.75,1.75];
            nurbs2{1}.coeffs(1:2,end,end) = nurbs2{2}.coeffs(1:2,1,end);
        case 5 % Two boxes with irregularities at one patch to be optimized (including a shared CP)
            nurbs = [nurbs,nurbs];
            nurbs(2) = translateNURBS(nurbs(2),Lvec);
            nurbs2 = nurbs;
            nurbs2{2}.coeffs(1:2,end,end) = [L+s,s];
            nurbs2{2}.coeffs(1:2,end-1,end) = [1.75,1.75];
            nurbs2{1}.coeffs(1:2,end,end) = nurbs2{2}.coeffs(1:2,end-1,end);
    end
end
% varCol.BC = 'None';    % not needed
% varCol.media = 'solid';
% varCol.analyticSolutionExist = false;
% varCol.applyLoad = 'None'; % no inhomogeneous Neumann condition
varCol.boundaryMethod = false;

if study > 5
    geometry = getTopology(nurbs2);
    topologysets = geometry.topologysets;
    [innerSurface,surf2volMap] = extractFreeSurface(nurbs2,'topologysets',topologysets,'topoSetName','inner');
    [outerSurface,outerSurf2volMap] = extractFreeSurface(nurbs2,'topologysets',topologysets,'topoSetName','outer');
end

% varCol.extraGP = 0;
% varCol.force = @(v,stopForDebug,J) -maxCompression(v,J,stopForDebug);
% varCol.C = elasticityMatrix(E,nu);
varCol.dimension = d;
varCol.nurbs = nurbs2;
% varCol.operator = 'linearElasticity';
varCol.fieldDimension = d;
varCol.extraGP = extraGP;
% varCol.applyBodyLoading = true;
% varCol.buildMassMatrix = false;
% varCol.progressBars = 0;
% varCol.buildStiffnessMatrix = true;

varCol = convertNURBS(varCol);
varCol = generateIGAmesh(varCol);
varCol = findDofsToRemove(varCol);
if study > 5
    varCol.geometry = geometry;
    varColBdry = meshBoundary(varCol,'inner');
    varColOuter = meshBoundary(varCol,'outer');
    nodes = varColOuter.nodes;
end
controlPts = varCol.controlPts;
x_cm = sum(controlPts,1)/size(controlPts,1);

noPatches = numel(nurbs2);
if study > 5
%     fixedBdryNodes = varColBdry.nodes; 
    fixedBdryNodes = [];
    fixedXdofs = varColBdry.nodes;
    fixedYdofs = varColBdry.nodes;
    fixedZdofs = varColBdry.nodes;
else
    if d == 2
        fixedXdofs = unique([find(abs(controlPts(:,1)) < Eps).', ...
                             find(abs(controlPts(:,2)) < Eps).']);
    else
        fixedXdofs = unique([find(abs(controlPts(:,1)) < Eps).', ...
                             find(abs(controlPts(:,2)) < Eps).', ...
                             find(abs(controlPts(:,3)) < Eps).']);     
        fixedZdofs = fixedXdofs;  
    end
    fixedYdofs = fixedXdofs;    
    fixedBdryNodes = [];
end
d = nurbs{1}.d;
if d == 2
    fixeddofs = union(fixedXdofs*d-1,fixedYdofs*d); 
else
    fixeddofs = union(union(fixedXdofs*d-2,fixedYdofs*d-1),fixedZdofs*d); 
end
varCol.dofsToRemove = union(union(fixeddofs,varCol.dofsToRemove),fixedBdryNodes);
varCol.freeDofs = setdiff(1:varCol.noDofs,varCol.dofsToRemove);
controlPtsT = controlPts.';

if prePlot    
    figure
    plotNURBS(nurbs2,'plotControlPolygon',true,'plotParmDir',0,'useLogForColoring', true,'plotAffinity',plotAffinity,'plotMeshQuality', plotMeshQuality, 'plotOrthogonality',plotOrthogonality,'resolution',resolution,'coarseLinearSampling',false)
    plotControlPts2(varCol,{fixedXdofs,  fixedYdofs,  fixedZdofs},...
                           {'fixedXdofs','fixedYdofs','fixedZdofs'})
    if d_p == 2
        view([0,90])
    else
        view([150,36])
%         view([77,17])
    end
    axis equal
    if d > 2
        camlight
    end
    colorbar
%     clim([0,1])
end

fig = figure;
switch method 
    case 'Nelder-Mead'
        options = optimset('MaxFunEvals', maxIter, 'MaxIter', maxIter, 'TolX',TolX,'TolFun',TolFun,'Display','iter');
        q = fminsearch(@(q) objFun(varCol, q, m, objectiveFun, useJacobian), controlPtsT(varCol.freeDofs), options);
    case 'Gradient decent'
        q = controlPtsT(varCol.freeDofs).';
        relativeDeltaq = Inf;
        gamma = 0.1;
        fprintf('\n%10s %12s %12s %12s\n', 'Iteration', 'f(q)', 'Rel. delta q', 'gamma');
        itercount = 0;
        [f, nabla_f] = objFun(varCol, q, m, objectiveFun, useJacobian);
        while relativeDeltaq > TolX && itercount < maxIter % && f > TolFun
            itercount = itercount + 1;
            fprintf('%10d %12.6g %12.6g %12.6g\n', itercount, f, relativeDeltaq, gamma);
            q_prev = q;
            nabla_f_prev = nabla_f;
    
            q = q - gamma*nabla_f;
    
            [f, nabla_f] = objFun(varCol, q, m, objectiveFun, useJacobian);
            gamma = abs((q-q_prev).'*(nabla_f - nabla_f_prev))/norm(nabla_f - nabla_f_prev).^2;
            relativeDeltaq = norm(q - q_prev)/norm(q);
            if plotIterations
                nurbs3 = updateNURBS(varCol,nurbs2,q);
                s = sprintf('Iteration = %3d\n', itercount);
                title(s)
                clf(fig)
                [~,minC,maxC] = plotNURBS(nurbs3,'plotControlPolygon',true,'plotParmDir',0,'useLogForColoring', true,'plotMeshQuality', plotMeshQuality, 'plotAffinity',plotAffinity,'plotOrthogonality',plotOrthogonality,'resolution',resolution,'coarseLinearSampling',false);
                if d_p == 2
                    view([0,90])
                else
                    view([77,17])
                end
                axis equal
                if d > 2
                    camlight
                end
                colorbar
                clim([0,1])
                s = sprintf('S2V/iterations/study%d/iteration%d.png', study, itercount);
                drawnow
%                 export_fig(s, '-png', '-r200')
    %             export_fig(s, '-png', '-transparent', '-r100')
                
            end
        end
end
nurbs3 = updateNURBS(varCol,nurbs2,q);

if postPlot
    figure
    title('Result after mesh optimization')
    plotNURBS(nurbs3,'plotControlPolygon',true,'plotParmDir',0,'useLogForColoring', true, 'plotMeshQuality', plotMeshQuality, 'plotAffinity',plotAffinity,'plotOrthogonality',plotOrthogonality,'resolution',resolution,'coarseLinearSampling',false)
    if d_p == 2
        view([0,90])
    else
        view([77,17])
    end
    axis equal
    colorbar
    if d > 2
        camlight
    end
    % clim([0,1])
%     clim([-17,-12])
end

if study < 6
    l2_error = 0;
    for patch = 1:noPatches
        l2_error = l2_error + sum(vecnorm(nurbs{patch}.coeffs(1:d,:) - nurbs3{patch}.coeffs(1:d,:), 2, 1).^2);
    end
    l2_error = sqrt(l2_error);
    fprintf('\nl2-error = %g\n', l2_error)
end

function nurbs = updateNURBS(varCol,nurbs,q)
controlPts = varCol.controlPts;
noPatches = numel(nurbs);
controlPtsT = controlPts.';
controlPtsT(varCol.freeDofs) = q;
gluedNodes = varCol.gluedNodes;
for i = 1:length(gluedNodes)
    parentIdx = gluedNodes{i}(1);
    for j = 2:length(gluedNodes{i})
        controlPtsT(:,gluedNodes{i}(j)) = controlPtsT(:,parentIdx);
    end
end
noPrevCP = 0;
d = nurbs{1}.d;
for patch = 1:noPatches
    indices_patch = (1:varCol.noCtrlPtsPatch(patch))+noPrevCP;
    nurbs{patch}.coeffs(1:d,:) = controlPtsT(:,indices_patch);
    noPrevCP = noPrevCP + varCol.noCtrlPtsPatch(patch);
end
end

function [objFunSum, nabla_objFunSum] = objFun(varCol, free_coeffs, m, objectiveFun, useJacobian)

compute_nabla_objFun = nargout > 1;
degree = varCol.degree;
knotVecs = varCol.knotVecs;
index = varCol.index;
noElems = varCol.noElems;
elRange = varCol.elRange;
element = varCol.element;
element2 = varCol.element2;
weights = varCol.weights;
controlPts = varCol.controlPts;
pIndex = varCol.pIndex;
noDofs = varCol.noDofs;
extraGP = varCol.extraGP;

controlPts = controlPts.';
controlPts(varCol.freeDofs) = free_coeffs.';
controlPts = controlPts.';

d = varCol.patches{1}.nurbs.d;
d_p = varCol.patches{1}.nurbs.d_p;
n_en = prod(degree+1);
[Q, W] = gaussTensorQuad(degree+1+extraGP(1:d_p));
% [Q, W] = gaussTensorQuad([1,1]);
noxi = numel(W);

objFunSum = 0; 
if compute_nabla_objFun
    nabla_objValues = zeros(d*n_en,noElems); 
    nabla_objIndices = zeros(d*n_en,noElems); 
    n = 2;
else
    if objectiveFun > 2
        n = 2;
    else
        n = 1;
    end
end
computeMixedDerivs = objectiveFun > 2;

% for e = 1:noElems
parfor e = 1:noElems
    patch = pIndex(e);
    knots = knotVecs{patch};
    Xi_e = zeros(d_p,2);
    for i = 1:d_p
        Xi_e(i,:) = elRange{i}(index(e,i),:);
    end

    sctr = element(e,:);
    pts = controlPts(sctr,:);
    wgts = weights(element2(e,:),:);
    
    J_2 = prod(Xi_e(:,2)-Xi_e(:,1))/2^d_p;
    
    sctr_k_e = zeros(1,d*n_en);
    for i = 1:d
        sctr_k_e(i:d:end) = d*(sctr-1)+i;
    end
    xi = parent2ParametricSpace(Xi_e, Q);
    I = findKnotSpans(degree, xi(1,:), knots);
    R = NURBSbasis(I, xi, degree, knots, wgts, n, computeMixedDerivs);

    fact = J_2.* W;
    if d_p == 2
        dXdxi = R{2}(:,:,1)*pts;
        dXdeta = R{3}(:,:,1)*pts;
        if computeMixedDerivs || compute_nabla_objFun
            d2Xdxi2 = R{2}(:,:,2)*pts;
            d2Xdeta2 = R{3}(:,:,2)*pts;
            if computeMixedDerivs
                d2Xdxieta = R{4}*pts;
            end
        end
        if useJacobian
            J_1 = dXdxi(:,1).*dXdeta(:,2) - dXdeta(:,1).*dXdxi(:,2);
        else
            J_1 = dot(dXdxi,dXdeta,2);
        end
    else
        dXdxi = R{2}(:,:,1)*pts;
        dXdeta = R{3}(:,:,1)*pts;
        dXdzeta = R{4}(:,:,1)*pts;
        if computeMixedDerivs || compute_nabla_objFun
            d2Xdxi2 = R{2}(:,:,2)*pts;
            d2Xdeta2 = R{3}(:,:,2)*pts;
            d2Xdzeta2 = R{4}(:,:,2)*pts;
        end
        if computeMixedDerivs
            d2Xdetazeta = R{5}*pts;
            d2Xdxizeta = R{6}*pts;
            d2Xdxieta = R{7}*pts;
        end
        if useJacobian
            J_1 = dot(dXdxi,cross(dXdeta,dXdzeta,2),2);
        else
            J_1 = dot(dXdeta,dXdzeta,2).^2+dot(dXdxi,dXdzeta,2).^2+dot(dXdxi,dXdeta,2).^2;
        end
    end

    if d_p == 2
        switch objectiveFun
            case 1
                sumdX = norm2(dXdxi).^d_p + norm2(dXdeta).^d_p;
                J_r = d_p*J_1./sumdX;
                if useJacobian
                    objFunSum = objFunSum + sum((J_r-1).^m.*fact,1);
                else
                    objFunSum = objFunSum + sum(J_r.^m.*fact,1);
                end
            case 2
                proddX = norm2(dXdxi).*norm2(dXdeta);
                J_r = J_1./proddX;
                if useJacobian
                    objFunSum = objFunSum + sum((J_r-1).^m.*fact,1);
                else
                    objFunSum = objFunSum + sum(J_r.^m.*fact,1);
                end
            case 3
                proddX = norm2(dXdxi).*norm2(dXdeta);
                sumd2X = dot(d2Xdxi2,d2Xdxi2,2) + dot(d2Xdeta2,d2Xdeta2,2) + 2*dot(d2Xdxieta,d2Xdxieta,2);
                J_r = J_1./proddX;
                if useJacobian
                    objFunSum = objFunSum + sum((J_r-1).^m.*sumd2X.*fact,1);
                else
                    objFunSum = objFunSum + sum(J_r.^m.*sumd2X.*fact,1);
                end
            case 4
                proddX = norm2(dXdxi).*norm2(dXdeta);
                sumd2X = dot(d2Xdxieta,d2Xdxieta,2);
                J_r = J_1./proddX;
                if useJacobian
                    objFunSum = objFunSum + sum((J_r-1).^m.*sumd2X.*fact,1);
                else
                    objFunSum = objFunSum + sum(J_r.^m.*sumd2X.*fact,1);
                end
            case 5
                sumd2X = dot(d2Xdxieta,d2Xdxieta,2);
                objFunSum = objFunSum + sum(2*sumd2X.*fact,1);
        end
        if compute_nabla_objFun
            d2Xdxidq = zeros(noxi,d,d*n_en);
            d2Xdetadq = zeros(noxi,d,d*n_en);
    
            d3Xdxi2dq = zeros(noxi,d,d*n_en);
            d3Xdeta2dq = zeros(noxi,d,d*n_en);
    
            d2Xdxidq(:,1,1:d:d*n_en) = R{2}(:,:,1);
            d2Xdxidq(:,2,2:d:d*n_en) = R{2}(:,:,1);
    
            d2Xdetadq(:,1,1:d:d*n_en) = R{3}(:,:,1);
            d2Xdetadq(:,2,2:d:d*n_en) = R{3}(:,:,1);
    
            d3Xdxi2dq(:,1,1:d:d*n_en) = R{2}(:,:,2);
            d3Xdxi2dq(:,2,2:d:d*n_en) = R{2}(:,:,2);
    
            d3Xdeta2dq(:,1,1:d:d*n_en) = R{3}(:,:,2);
            d3Xdeta2dq(:,2,2:d:d*n_en) = R{3}(:,:,2);

            if computeMixedDerivs
                d3Xdxietadq = zeros(noxi,d,d*n_en);
                d3Xdxietadq(:,1,1:d:d*n_en) = R{4};
                d3Xdxietadq(:,2,2:d:d*n_en) = R{4};
    
                d2XdxietaRep = repmat(d2Xdxieta,1,1,d*n_en);
            end
    
            dXdxiRep = repmat(dXdxi,1,1,d*n_en);
            dXdetaRep = repmat(dXdeta,1,1,d*n_en);
    
            d2Xdxi2Rep = repmat(d2Xdxi2,1,1,d*n_en);
            d2Xdeta2Rep = repmat(d2Xdeta2,1,1,d*n_en);
    
            if useJacobian
                nabla_J_1 = d2Xdxidq(:,1,:).*dXdeta(:,2) + dXdxi(:,1).*d2Xdetadq(:,2,:) - (d2Xdetadq(:,1,:).*dXdxi(:,2) + dXdeta(:,1).*d2Xdxidq(:,2,:));
            else
                nabla_J_1 = dot(d2Xdxidq,dXdetaRep,2) + dot(dXdxiRep,d2Xdetadq,2);
            end
            switch objectiveFun
                case 1
                    sumdX2 = norm2(dXdxi).^(d_p-2).*dot(d2Xdxidq,dXdxiRep,2) + norm2(dXdeta).^(d_p-2).*dot(d2Xdetadq,dXdetaRep,2);
                    nabla_J_r = d_p*nabla_J_1./sumdX - d_p^2*J_1./sumdX.^2.*sumdX2;
                    if useJacobian
                        nabla_objValues(:,e) = sum(m.*(J_r-1).^(m-1).*nabla_J_r.*fact, 1);
                    else
                        nabla_objValues(:,e) = sum(m.*J_r.^(m-1).*nabla_J_r.*fact, 1);
                    end
                case 2
                    sumdX2 =  norm2(dXdeta).*dot(d2Xdxidq,dXdxiRep,2)./norm2(dXdxi) ...
                            + norm2(dXdxi).*dot(d2Xdetadq,dXdetaRep,2)./norm2(dXdeta);
                    nabla_J_r = nabla_J_1./proddX - J_1./proddX.^2.*sumdX2;
                    if useJacobian
                        nabla_objValues(:,e) = sum(m.*(J_r-1).^(m-1).*nabla_J_r.*fact, 1);
                    else
                        nabla_objValues(:,e) = sum(m.*J_r.^(m-1).*nabla_J_r.*fact, 1);
                    end
                case 3
                    sumdX2 =  norm2(dXdeta).*dot(d2Xdxidq,dXdxiRep,2)./norm2(dXdxi) ...
                            + norm2(dXdxi).*dot(d2Xdetadq,dXdetaRep,2)./norm2(dXdeta);
                    sumd2X2 = dot(d3Xdxi2dq,d2Xdxi2Rep,2) + dot(d3Xdeta2dq,d2Xdeta2Rep,2) + 2*dot(d3Xdxietadq,d2XdxietaRep,2);
                    nabla_J_r = nabla_J_1./proddX - J_1./proddX.^2.*sumdX2;
                    if useJacobian
                        nabla_objValues(:,e) = sum((m.*(J_r-1).^(m-1).*nabla_J_r.*sumd2X + 2*(J_r-1).^m.*sumd2X2).*fact, 1);
                    else
                        nabla_objValues(:,e) = sum((m.*J_r.^(m-1).*nabla_J_r.*sumd2X + 2*J_r.^m.*sumd2X2).*fact, 1);
                    end
                case 4
                    sumdX2 =  norm2(dXdeta).*dot(d2Xdxidq,dXdxiRep,2)./norm2(dXdxi) ...
                            + norm2(dXdxi).*dot(d2Xdetadq,dXdetaRep,2)./norm2(dXdeta);
                    sumd2X2 = dot(d3Xdxietadq,d2XdxietaRep,2);
                    nabla_J_r = nabla_J_1./proddX - J_1./proddX.^2.*sumdX2;
                    if useJacobian
                        nabla_objValues(:,e) = sum((m.*(J_r-1).^(m-1).*nabla_J_r.*sumd2X + 2*(J_r-1).^m.*sumd2X2).*fact, 1);
                    else
                        nabla_objValues(:,e) = sum((m.*J_r.^(m-1).*nabla_J_r.*sumd2X + 2*J_r.^m.*sumd2X2).*fact, 1);
                    end
                case 5
                    sumd2X2 = dot(d3Xdxietadq,d2XdxietaRep,2);
                    nabla_objValues(:,e) = sum(sumd2X2.*fact, 1);
            end
    
            nabla_objIndices(:,e) = sctr_k_e;
        end
    else
        if ~useJacobian
            error('Not implemented')
        end
        flag = 0;
        switch objectiveFun
            case 1
                sumdX = norm2(dXdxi).^d_p + norm2(dXdeta).^d_p + norm2(dXdzeta).^d_p;
                J_r = d_p*J_1./sumdX;
                objFunSum = objFunSum + sum((J_r-1).^m.*fact,1);
            case 2
                proddX = norm2(dXdxi).*norm2(dXdeta).*norm2(dXdzeta);
                J_r = J_1./proddX;
                objFunSum = objFunSum + sum((J_r-1).^m.*fact,1);
            case 3
                proddX = norm2(dXdxi).*norm2(dXdeta).*norm2(dXdzeta);
                sumd2X =   dot(d2Xdxi2,d2Xdxi2,2) + dot(d2Xdeta2,d2Xdeta2,2) + dot(d2Xdzeta2,d2Xdzeta2,2) ...
                         + 2*dot(d2Xdetazeta,d2Xdetazeta,2) + 2*dot(d2Xdxizeta,d2Xdxizeta,2) + 2*dot(d2Xdxieta,d2Xdxieta,2);
                if flag
                    J_r = J_1./proddX.*sumd2X;
    %                 objFunSum = objFunSum + sum((J_r-1).^m.*fact,1);
                    objFunSum = objFunSum + sum((J_1./proddX-1).^m.*sumd2X.*fact,1);
                else
                    J_r = J_1./proddX;
                    objFunSum = objFunSum + sum((J_r-1).^m.*sumd2X.*fact,1);
                end
            case 4
                proddX = norm2(dXdxi).*norm2(dXdeta).*norm2(dXdzeta);
                sumd2X = dot(d2Xdetazeta,d2Xdetazeta,2) + dot(d2Xdxizeta,d2Xdxizeta,2) + dot(d2Xdxieta,d2Xdxieta,2);
                J_r = J_1./proddX;
                objFunSum = objFunSum + sum((J_r-1).^m.*sumd2X.*fact,1);
            case 5
                sumd2X = dot(d2Xdetazeta,d2Xdetazeta,2) + dot(d2Xdxizeta,d2Xdxizeta,2) + dot(d2Xdxieta,d2Xdxieta,2);
                objFunSum = objFunSum + sum(sumd2X.*fact,1);
        end
        if compute_nabla_objFun
            d2Xdxidq = zeros(noxi,d,d*n_en);
            d2Xdetadq = zeros(noxi,d,d*n_en);
            d2Xdzetadq = zeros(noxi,d,d*n_en);
    
            d3Xdxi2dq = zeros(noxi,d,d*n_en);
            d3Xdeta2dq = zeros(noxi,d,d*n_en);
            d3Xdzeta2dq = zeros(noxi,d,d*n_en);
    
            d2Xdxidq(:,1,1:d:d*n_en) = R{2}(:,:,1);
            d2Xdxidq(:,2,2:d:d*n_en) = R{2}(:,:,1);
            d2Xdxidq(:,3,3:d:d*n_en) = R{2}(:,:,1);
    
            d2Xdetadq(:,1,1:d:d*n_en) = R{3}(:,:,1);
            d2Xdetadq(:,2,2:d:d*n_en) = R{3}(:,:,1);
            d2Xdetadq(:,3,3:d:d*n_en) = R{3}(:,:,1);
    
            d2Xdzetadq(:,1,1:d:d*n_en) = R{4}(:,:,1);
            d2Xdzetadq(:,2,2:d:d*n_en) = R{4}(:,:,1);
            d2Xdzetadq(:,3,3:d:d*n_en) = R{4}(:,:,1);
    
            d3Xdxi2dq(:,1,1:d:d*n_en) = R{2}(:,:,2);
            d3Xdxi2dq(:,2,2:d:d*n_en) = R{2}(:,:,2);
            d3Xdxi2dq(:,3,3:d:d*n_en) = R{2}(:,:,2);
    
            d3Xdeta2dq(:,1,1:d:d*n_en) = R{3}(:,:,2);
            d3Xdeta2dq(:,2,2:d:d*n_en) = R{3}(:,:,2);
            d3Xdeta2dq(:,3,3:d:d*n_en) = R{3}(:,:,2);
    
            d3Xdzeta2dq(:,1,1:d:d*n_en) = R{4}(:,:,2);
            d3Xdzeta2dq(:,2,2:d:d*n_en) = R{4}(:,:,2);
            d3Xdzeta2dq(:,3,3:d:d*n_en) = R{4}(:,:,2);
        
            dXdxiRep = repmat(dXdxi,1,1,d*n_en);
            dXdetaRep = repmat(dXdeta,1,1,d*n_en);
            dXdzetaRep = repmat(dXdzeta,1,1,d*n_en);
    
            d2Xdxi2Rep = repmat(d2Xdxi2,1,1,d*n_en);
            d2Xdeta2Rep = repmat(d2Xdeta2,1,1,d*n_en);
            d2Xdzeta2Rep = repmat(d2Xdzeta2,1,1,d*n_en);
    
            if computeMixedDerivs    
                d3Xdetazetadq = zeros(noxi,d,d*n_en);
                d3Xdxizetadq = zeros(noxi,d,d*n_en);
                d3Xdxietadq = zeros(noxi,d,d*n_en);
                
                d3Xdetazetadq(:,1,1:d:d*n_en) = R{5};
                d3Xdetazetadq(:,2,2:d:d*n_en) = R{5};
                d3Xdetazetadq(:,3,3:d:d*n_en) = R{5};
        
                d3Xdxizetadq(:,1,1:d:d*n_en) = R{6};
                d3Xdxizetadq(:,2,2:d:d*n_en) = R{6};
                d3Xdxizetadq(:,3,3:d:d*n_en) = R{6};
        
                d3Xdxietadq(:,1,1:d:d*n_en) = R{7};
                d3Xdxietadq(:,2,2:d:d*n_en) = R{7};
                d3Xdxietadq(:,3,3:d:d*n_en) = R{7};
    
                d2XdetazetaRep = repmat(d2Xdetazeta,1,1,d*n_en);
                d2XdxizetaRep = repmat(d2Xdxizeta,1,1,d*n_en);
                d2XdxietaRep = repmat(d2Xdxieta,1,1,d*n_en);
            end
    
            nabla_J_1 = dot(d2Xdxidq,cross(dXdetaRep,dXdzetaRep,2),2) + dot(dXdxiRep,cross(d2Xdetadq,dXdzetaRep,2) + cross(dXdetaRep,d2Xdzetadq,2),2);
            switch objectiveFun
                case 1
                    sumdX2 = norm2(dXdxi).^(d_p-2).*dot(d2Xdxidq,dXdxiRep,2) + norm2(dXdeta).^(d_p-2).*dot(d2Xdetadq,dXdetaRep,2) + norm2(dXdzeta).^(d_p-2).*dot(d2Xdzetadq,dXdzetaRep,2);
                    nabla_J_r = d_p*nabla_J_1./sumdX - d_p^2*J_1./sumdX.^2.*sumdX2;
                    nabla_objValues(:,e) = sum(m.*(J_r-1).^(m-1).*nabla_J_r.*fact, 1);
                case 2
                    sumdX2 =  norm2(dXdeta).*norm2(dXdzeta).*dot(d2Xdxidq,dXdxiRep,2)./norm2(dXdxi) ...
                            + norm2(dXdxi).*norm2(dXdzeta).*dot(d2Xdetadq,dXdetaRep,2)./norm2(dXdeta) ...
                            + norm2(dXdxi).*norm2(dXdeta).*dot(d2Xdzetadq,dXdzetaRep,2)./norm2(dXdzeta);
                    nabla_J_r = nabla_J_1./proddX - J_1./proddX.^2.*sumdX2;
                    nabla_objValues(:,e) = sum(m.*(J_r-1).^(m-1).*nabla_J_r.*fact, 1);
                case 3
                    sumdX2 =  norm2(dXdeta).*norm2(dXdzeta).*dot(d2Xdxidq,dXdxiRep,2)./norm2(dXdxi) ...
                            + norm2(dXdxi).*norm2(dXdzeta).*dot(d2Xdetadq,dXdetaRep,2)./norm2(dXdeta) ...
                            + norm2(dXdxi).*norm2(dXdeta).*dot(d2Xdzetadq,dXdzetaRep,2)./norm2(dXdzeta);
    %                 sumd2X =   dot(d2Xdxi2,d2Xdxi2,2) + dot(d2Xdeta2,d2Xdeta2,2) + dot(d2Xdzeta2,d2Xdzeta2,2) ...
    %                          + 2*dot(d2Xdetazeta,d2Xdetazeta,2) + 2*dot(d2Xdxizeta,d2Xdxizeta,2) + 2*dot(d2Xdxieta,d2Xdxieta,2);
                    sumd2X2 =   dot(d3Xdxi2dq,d2Xdxi2Rep,2) + dot(d3Xdeta2dq,d2Xdeta2Rep,2) + dot(d3Xdzeta2dq,d2Xdzeta2Rep,2) + ...
                              + 2*dot(d3Xdetazetadq,d2XdetazetaRep,2) + 2*dot(d3Xdxizetadq,d2XdxizetaRep,2) + 2*dot(d3Xdxietadq,d2XdxietaRep,2);
                    nabla_J_r = nabla_J_1./proddX - J_1./proddX.^2.*sumdX2;
                    if flag
                        J_r = J_1./proddX.*sumd2X;
                        objFunSum = objFunSum + sum((J_r-1).^m.*fact,1);
                        nabla_objValues(:,e) = sum((m.*(J_r-1).^(m-1).*(nabla_J_r.*sumd2X + 2*J_1./proddX.*sumd2X2)).*fact, 1);
    %                     nabla_objValues(:,e) = sum((m.*(J_r-1).^(m-1).*nabla_J_r + 2*(J_r-1).^m.*sumd2X2).*fact, 1);
                    else
                        nabla_objValues(:,e) = sum((m.*(J_r-1).^(m-1).*nabla_J_r.*sumd2X + 2*(J_r-1).^m.*sumd2X2).*fact, 1);
                    end
                case 4
                    sumdX2 =  norm2(dXdeta).*norm2(dXdzeta).*dot(d2Xdxidq,dXdxiRep,2)./norm2(dXdxi) ...
                            + norm2(dXdxi).*norm2(dXdzeta).*dot(d2Xdetadq,dXdetaRep,2)./norm2(dXdeta) ...
                            + norm2(dXdxi).*norm2(dXdeta).*dot(d2Xdzetadq,dXdzetaRep,2)./norm2(dXdzeta);
                    sumd2X2 = dot(d3Xdetazetadq,d2XdetazetaRep,2) + dot(d3Xdxizetadq,d2XdxizetaRep,2) + dot(d3Xdxietadq,d2XdxietaRep,2);
                    nabla_J_r = nabla_J_1./proddX - J_1./proddX.^2.*sumdX2;
                    nabla_objValues(:,e) = sum((m.*(J_r-1).^(m-1).*nabla_J_r.*sumd2X + 2*(J_r-1).^m.*sumd2X2).*fact, 1);
                case 5
                    sumd2X2 = dot(d3Xdetazetadq,d2XdetazetaRep,2) + dot(d3Xdxizetadq,d2XdxizetaRep,2) + dot(d3Xdxietadq,d2XdxietaRep,2);
                    nabla_objValues(:,e) = sum(2*sumd2X2.*fact, 1);
            end
    
            nabla_objIndices(:,e) = sctr_k_e;
        end
    end
end
if compute_nabla_objFun
    nabla_objFunSum = vectorAssembly(nabla_objValues,nabla_objIndices,noDofs);
    nabla_objFunSum(varCol.dofsToRemove) = [];
    nabla_objFunSum = nabla_objFunSum(:);
end
end
%     
% function [objFunSum, nabla_objFunSum] = objFun(nurbs, noPatches, XIETAZETA, indices_free, free_coeffs)
% 
% compute_nabla_objFun = nargout > 1;
% if compute_nabla_objFun
%     nabla_objFunSum = zeros(1,1,numel(free_coeffs)); 
%     noxi = size(XIETAZETA,1);
% end
% % nurbs = updateGluedNodes()
% objFunSum = 0; 
% for patch = 1:noPatches
%     nurbs_patch = nurbs{patch};
%     indices_free_patch = indices_free{patch,1};
%     indices_free_glob = indices_free{patch,2};
%     nurbs_patch.coeffs(1:3,indices_free_patch) = free_coeffs(:,indices_free_glob);
% 
%     if compute_nabla_objFun
%         d_p = nurbs_patch.d_p;
%         d = nurbs_patch.d;
%         noDofs_patch = d*numel(indices_free_patch);
%         [~,dXdxi,dXdeta,dXdzeta,RdR,w_i] = evaluateNURBS(nurbs_patch,XIETAZETA,1);
%         J_1 = dot(dXdxi,cross(dXdeta,dXdzeta,2),2);
%         d2Xdxidq = zeros(noxi,d,noDofs_patch);
%         d2Xdetadq = zeros(noxi,d,noDofs_patch);
%         d2Xdzetadq = zeros(noxi,d,noDofs_patch);
% 
%         w_i = reshape(w_i,noxi,[]);
%         for j = 1:noDofs_patch
%             j_cp = ceil(j/d);
%             [val,j_loc] = max(indices_free_patch(j_cp) == w_i, [],2);
%             if val ~= 0
%                 i_e = mod(j-1,d)+1;
%                 indices = (j_loc-1)*noxi+(1:noxi).';
%                 d2Xdxidq(:,i_e,j) = RdR{2}(indices);
%                 d2Xdetadq(:,i_e,j) = RdR{3}(indices);
%                 d2Xdzetadq(:,i_e,j) = RdR{4}(indices);
%             end
%         end
% 
%         dXdxiRep = repmat(dXdxi,1,1,noDofs_patch);
%         dXdetaRep = repmat(dXdeta,1,1,noDofs_patch);
%         dXdzetaRep = repmat(dXdzeta,1,1,noDofs_patch);
%         sumdX = norm2(dXdxi).^d_p + norm2(dXdeta).^d_p + norm2(dXdzeta).^d_p;
%         sumdX2 = norm2(dXdxi).^(d_p-2).*dot(d2Xdxidq,dXdxiRep,2) + norm2(dXdeta).^(d_p-2).*dot(d2Xdetadq,dXdetaRep,2) + norm2(dXdzeta).^(d_p-2).*dot(d2Xdzetadq,dXdzetaRep,2);
%         J_r = d_p*J_1./sumdX;
%         nabla_J_1 = dot(d2Xdxidq,cross(dXdetaRep,dXdzetaRep,2),2) + dot(dXdxiRep,cross(d2Xdetadq,dXdzetaRep,2) + cross(dXdetaRep,d2Xdzetadq,2),2);
%         nabla_J_r = d_p*nabla_J_1./sumdX - d_p^2*J_1./sumdX.^2.*sumdX2;
%     
%         indices_free_glob2 = [indices_free_glob*d-2; indices_free_glob*d-1; indices_free_glob*d];
%         nabla_objFunSum(1,1,indices_free_glob2(:)) = nabla_objFunSum(1,1,indices_free_glob2(:)) + 2*sum((J_r-1).*nabla_J_r, 1);
%     else
%         J_r = meanRatioJacobian(nurbs_patch,XIETAZETA);
%     end
%     objFunSum = objFunSum + norm(J_r-1).^2;
% end
% objFunSum = objFunSum/noPatches;
% if compute_nabla_objFun
%     nabla_objFunSum = nabla_objFunSum(:)/noPatches;
% end
% end
% addpath(genpath('../export_fig'))
% no2Dpoints = 1000;
% startMatlabPool
% para.plotResultsInParaview = 1;
% clear varCol
% model = 'cube';
% prePlot.plot3Dgeometry = 0;
% homeDir = expanduser('~');
% 
% 
% 
% for M = 1
%     switch model
%         case 'cube'
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %% Make Jacobian example
%             close all
%             nurbs = getPrismData();
%             nurbs = elevateNURBSdegree(nurbs,1);
%             nurbs2 = nurbs;
%             s = 0.7;
%             nurbs2{1}.coeffs([1,3],end,end-1,end) = [s,s];
%             nurbs2{1}.coeffs([2,3],end-1,end,end) = [s,s];
%             nurbs2{1}.coeffs([1,2],end,end,end-1) = [s,s];
%             
%             %% We want to make an mesh optimizer that recreates nurbs from nurbs2
% %             nurbs = insertKnotsInNURBS(nurbs,[1,2,3]);
% %             nurbs2 = insertKnotsInNURBS(nurbs2,[1,2,3]);
%             if 0
%                 plotNURBS(nurbs,'plotControlPolygon',true,'plotParmDir',0,'plotJacobian',true,'resolution',[100,100,100])
%                 view([77,17])
%                 axis equal
%                 camlight
%                 colorbar
%                 clim([0,1])
%                 
%                 figure
%                 plotNURBS(nurbs2,'plotControlPolygon',true,'plotParmDir',0,'plotJacobian',true,'resolution',[100,100,100])
%                 view([77,17])
%                 axis equal
%                 camlight
%                 colorbar
%                 clim([0,1])
%             end
%             nurbs = nurbs2;
% 
%             E = 3.0e10;
%             nu = 0.2;
%             rho = 2500;
%             noInsPts = 3;
%             newKnots = (2^(M-1)-1)*[1,1,1];
%             paraExtraScale = 20;
%     end
%     varCol.BC = 'None';    % not needed
%     varCol.media = 'solid';
%     varCol.analyticSolutionExist = false;
%     varCol.applyLoad = 'None'; % no inhomogeneous Neumann condition
%     varCol.boundaryMethod = false;
%     nurbs = insertKnotsInNURBS(nurbs,newKnots);
% 
%     varCol.extraGP = 0;
%     varCol.force = @(v,stopForDebug,J) maxCompression(v,J,stopForDebug);
%     varCol.C = elasticityMatrix(E,nu);
%     varCol.dimension = 3;
%     varCol.nurbs = nurbs;
%     varCol.operator = 'linearElasticity';
%     varCol.fieldDimension = 3;
%     varCol.applyBodyLoading = true;
%     varCol.buildMassMatrix = false;
%     varCol.progressBars = 0;
%     varCol.buildStiffnessMatrix = true;
%     
%     varCol = convertNURBS(varCol);
%     varCol = generateIGAmesh(varCol);
%     varCol = findDofsToRemove(varCol);
%     noDofs = varCol.noDofs;
% 
%     stringShift = 40;
%     controlPts = varCol.controlPts;
%     Eps = 1e-7;
%     if prePlot.plot3Dgeometry
%         plotNURBS(nurbs,'color',getColor(1),'resolution',30*ones(1,3),'plotControlPolygon',1,'plotJacobian',true)
%         axis equal
%         view(getView())
%         camlight
%         colorbar
%         clim([0,1])
%     end
%     switch model
%         case 'cube'
%             fixedXdofs = find(abs(controlPts(:,1)) < Eps).';        
%             fixedYdofs = find(abs(controlPts(:,2)) < Eps).';        
%             fixedZdofs = find(abs(controlPts(:,3)) < Eps).';    
%             rigiddofs = unique([fixedXdofs,fixedYdofs,fixedZdofs]);
%             if prePlot.plot3Dgeometry
%                 plotControlPts2(varCol,{rigiddofs,  rigiddofs,  rigiddofs},...
%                                        {'fixedXdofs','fixedYdofs','fixedZdofs'})
%             end
%     end
% 
%     varCol.progressBars = true;
%     d = varCol.dimension;
%     fixeddofs = union(union(rigiddofs*d-2,rigiddofs*d-1),rigiddofs*d); 
%     dofsToRemove = varCol.dofsToRemove;
%     dofsToRemove = union(fixeddofs,dofsToRemove);
% 
%     %% Build stiffness matrix
%     fprintf(['\n%-' num2str(stringShift) 's'], 'Building system matrices ... ')
%     tic
%     task.varCol = varCol;
%     task.misc.extraGP = [0,0,0];
%     task.analyticSolutionExist = false;
%     task.splitExteriorFields = false;
%     task.misc.BC = 'SHBC';
%     task = buildMatrices(task);
%     A = task.varCol.A_K;
%     F = task.varCol.FF;
%     fprintf('using %12f seconds.', toc)
% 
%     %% Modify stiffness matrix and force vector due to glued nodes and homogeneous dirichlet conditions
%     A(dofsToRemove,:) = [];
%     A(:,dofsToRemove) = [];   
%     F(dofsToRemove) = [];  
%     fprintf(['\n%-' num2str(stringShift) 's'], 'Solving system of equations ... ') 
%     tic
%     UU = A\F;
%     fprintf('using %12f seconds.', toc)
%     U = zeros(noDofs,1);
%     U(setdiff(1:noDofs, dofsToRemove'),:) = UU;   
%     task.varCol.U = U;
%     task.varCol = addSolutionToRemovedNodes(task.varCol);
%     para.U = task.varCol.U;
%     para.name = [homeDir '/results/ASIGA/' model '/' model '_M' num2str(M)];
%     if para.plotResultsInParaview
% 
%         task.varCol.isOuterDomain = false;
%         d = task.varCol.dimension;
%         isSolid = d == 3;
% 
%         para.celltype = 'VTK_HEXAHEDRON';
%         para.plotP_inc = false;
%         para.plotScalarField = false;
%         para.plotTotField = false; 
%         para.plotTotFieldAbs = false; 
%         para.plotAnalytic = false; 
%         para.plotTimeOscillation = false;
%         para.computeGrad = false;
%         para.plotVonMisesStress = true;
%         para.plotStressXX = true;
%         para.plotStressYY = true;
%         para.plotStressZZ = true;
%         para.plotStressYZ = true;
%         para.plotStressXZ = true;
%         para.plotStressXY = true;
%         para.plotError = false; 
%         para.plotArtificialBoundary = false;
%         para.plotDisplacementVectors = true;
% 
%         para.extraXiPts   = round(paraExtraScale/2^(M-1));  % Extra visualization points in the xi-direction per element
%         para.extraEtaPts  = round(paraExtraScale/2^(M-1));  % Extra visualization points in the eta-direction per element
%         para.extraZetaPts = round(paraExtraScale/2^(M-1));  % Extra visualization points in the zeta-direction per element
% 
%         para.i_MS = 1;
%         para.plotMesh = true;
%         fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting results in Paraview ... ')
%         createParaviewFiles(task, 'U', {U}, 'para_options', para);
%         fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting mesh files in Paraview ... ')
%         tic
%         createVTKmeshFiles(task, 'U', {U}, 'para_options', para)
%         fprintf('using %12f seconds.\n', toc)
%     end
% end

    
function plotControlPts2(varCol,dofs,labels)

controlPts = varCol.controlPts;
noTypes = numel(dofs);
cmap = [1,0,0;
        1,1,0;
        0,0,1;
        1,0,1;
        0,0,0.5;
        0,1,1;
        0,1,0];

hold on
for i = 1:noTypes
    indices = dofs{i};
    p(i) = plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o',...
                'color',cmap(i,:),'MarkerFaceColor',cmap(i,:), 'MarkerEdgeColor', cmap(i,:));
end
legend(p,labels)
% legend show
end

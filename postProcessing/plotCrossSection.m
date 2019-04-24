function plotCrossSection(varCol, U, options, parm_pt, dir, delta, xb, yb, zb, useExtraQuadPts, nameAppendix)


Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);

weights = varCol.weights;
controlPts = varCol.controlPts;

runInParallell = varCol.runInParallell;
[visElements, nodes, pts, pts2] = triangulateCrossSection(varCol.nurbs, parm_pt, dir, delta, xb, yb, zb);

if ~options.plotTimeOscillation
    options.plotTotFieldAbs = 1;
end

options.name = [options.name, nameAppendix];
options.celltype = 'VTK_TRIANGLE';

data.nodes = nodes;
data.visElements = visElements;
noKnots = length(pts);

totFieldOn = zeros(noKnots, 1);
v = zeros(noKnots,3);   
parfor (i = 1:noKnots, runInParallell)
% for i = 1:noKnots
    switch dir
        case 1
            [totFieldOn(i), v(i,:)] = numericalSolEval_final_surf(pts(i), parm_pt(2), p_xi, p_eta, Xi, Eta, weights, controlPts, U);
        case 2
            [totFieldOn(i), v(i,:)] = numericalSolEval_final_surf(parm_pt(1), pts(i), p_xi, p_eta, Xi, Eta, weights, controlPts, U);
    end    
end
if dir == 2
    noKnots2 = length(pts2);
    totFieldOn2 = zeros(noKnots2, 1);
    v2 = zeros(noKnots2,3);   
    parfor (i = 1:noKnots2, runInParallell)
    % for i = 1:noKnots
        [totFieldOn2(i), v2(i,:)] = numericalSolEval_final_surf(parm_pt(2), pts2(i), p_xi, p_eta, Xi, Eta, weights, controlPts, U);

    end
    totFieldOn = [totFieldOn; totFieldOn2];
    v = [v; v2];
    noKnots = noKnots + noKnots2;
end
p_inc = varCol.p_inc;
switch varCol.method
    case 'BEM'
        scalarField = calculateScatteredPressureBEM(varCol, U, nodes(noKnots+1:end,:), useExtraQuadPts, 0);
        scalarField = [totFieldOn-p_inc(v); scalarField];
        totField = scalarField + p_inc(nodes);
    case 'IE'
        scalarField = calculateScatteredPressure(varCol.varColFull, varCol.U_full, nodes, 0, 0);
        totField = scalarField + p_inc(nodes);
    case 'IENSG'
        scalarField = calculateScatteredPressureNonSepGeom(varCol.varColFull, varCol.U_full, nodes(noKnots+1:end,:), useExtraQuadPts, 0);
        scalarField = [totFieldOn; scalarField];
        totField = scalarField + p_inc(nodes);
end

omega = varCol.omega;
data.P_inc = real(makeDynamic(p_inc(nodes), options, omega)); 
data.scalarField = real(makeDynamic(scalarField, options, omega)); 
if options.plotAnalytic
    analytic = varCol.analytic(nodes);
    data.analytic = real(makeDynamic(analytic, options, omega)); 
    data.Error = abs(scalarField-analytic)./abs(analytic);
end
data.totField = real(makeDynamic(totField, options, omega));
data.totFieldAbs = abs(makeDynamic(totField, options, omega));

data.omega = omega;
options.plotDisplacementVectors = false;
makeVTKfile(data, options);







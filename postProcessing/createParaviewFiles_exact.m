function createParaviewFiles_exact(varCol, extraXiPts, extraEtaPts, extraZetaPts, options, P_inc, gP_inc, analytic, gAnalytic, analytic_solid, analytic_stress)

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
nurbs = varCol.nurbs;
% noElems = varCol.noElems;
d = varCol.dimension;

omega = varCol.omega;

if d == 1
%     [nodes, ~, visElements] = buildVisualization3dMesh_new3(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, nurbs);
    [nodes, ~, visElements] = buildVisualization3dCuttingPlaneMesh(Eta, Zeta, extraEtaPts, extraZetaPts, nurbs);
else
    [nodes, ~, visElements] = buildVisualization3dSurfaceMesh(Xi, Eta, extraXiPts, extraEtaPts, nurbs);
end
% [nodes, noNodes, visElements] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, nurbs);

data.nodes = nodes;
data.visElements = visElements;
if d == 1
    rho_f = varCol.rho_f;
    gScalarField = gAnalytic(nodes);
    displacement = (gScalarField+gP_inc(nodes))/(rho_f*omega^2);
    data.P_inc = P_inc(nodes);
    scalarField = analytic(nodes);
    data.totField = scalarField + P_inc(nodes);
    data.scalarField = scalarField;
    data.gradient = gScalarField;
else
    displacement = analytic_solid(nodes);
    data.stress = analytic_stress(nodes);
end
if max(max(abs(displacement))) < 1e-6
    data.displacement = 0.25*displacement*4e9;
else
    data.displacement = displacement;
end

data.omega = omega;
% data.noElems = noElems;
makeVTKfile(data, options);
% makeVTFfile(data, options);




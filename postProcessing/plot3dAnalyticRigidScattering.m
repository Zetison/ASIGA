

noUniqueXiKnots = length(unique(varCol_fluid_o.nurbs.knots{1}));
noUniqueEtaKnots = length(unique(varCol_fluid_o.nurbs.knots{2}));
noUniqueZetaKnots = length(unique(varCol_fluid_o.nurbs.knots{3}));

extraXiPts = floor(40/(noUniqueXiKnots-1)); % .. per element
extraEtaPts = floor(40/(noUniqueEtaKnots-1)); % .. per element
extraZetaPts = floor(0/(noUniqueZetaKnots-1)); % .. per element
                
vtfFileName = ['../graphics/paraview/' model '_' BC '/mesh' num2str(mesh) scatteringCase '_f_' num2str(freq) 'alpha' num2str(alpha_s_arr*180/pi) 'totalField0_outer_analytic_'];

options = {'name',vtfFileName, 'plotScalarField',1, 'plotTimeOscillation', 1};

varCol = varCol_fluid_o;

Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
Zeta = varCol.nurbs.knots{3};
p = varCol.nurbs.degree(1);
q = varCol.nurbs.degree(2);
r = varCol.nurbs.degree(3);

index = varCol.index;
noElems = varCol.noElems;
elRangeXi = varCol.elRangeXi;
elRangeEta = varCol.elRangeEta;
elRangeZeta = varCol.elRangeZeta;
element = varCol.element;
weights = varCol.weights;
controlPts = varCol.controlPts;
noCtrlPts = varCol.noCtrlPts;
noDofs = varCol.noDofs;
nurbs = varCol.nurbs;
if varCol.dimension == 3
    C = varCol.C;
    Ux = U(1:noCtrlPts);
    Uy = U(noCtrlPts+1:2*noCtrlPts);
    Uz = U(2*noCtrlPts+1:noDofs);
end
% build visualization B8 mesh

[nodes, noNodes, visElements, cornerNode, ...
          noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, nurbs);


displacement   = zeros(noNodes,3);

x = nodes(:,1);
y = nodes(:,2);
z = nodes(:,3);
% radius = sqrt(x.^2+y.^2+z.^2);
% theta = acos(z./radius);

% [r_m, theta_m, phi_m] = coordTransform(x,y,z,pi/4,pi/4);

scalarField = scatteredPressureOnRigidSphere2(x, y, z, P_0, k_wn_o, R_o, 1e-13,alpha_s_arr,beta_s_arr);


dataAnal.nodes = nodes;
dataAnal.visElements = visElements;
dataAnal.omega = omega;
dataAnal.scalarField = scalarField;

makeVTKfile_new(dataAnal, options);


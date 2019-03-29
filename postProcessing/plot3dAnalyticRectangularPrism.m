close all
clear all
for analyticFunction = 1:6
    analyticFunction
    if analyticFunction <= 2
        wx = 1;  %1
        wy = 0.01;  %0.01
        wz = 0.1; %0.1
    else
        wx = 3;
        wy = 2;
        wz = 1;
    end
    nu = 0.3;           % Poisson's ratio
    E = 13;           % Young modulus (in Pa)
    lambda = nu*E/(1+nu)/(1-2*nu);
    mu = E/(2*(1+nu));
    getAnalyticFunctions_RectangularPrism

    solid = getRectangularPrismData(wx, wy, wz, [wx/2, wy/2, wz/2]);

    [Xi, Eta, Zeta, ...
     p, q, r, ...
     n, m, l, ...
     weights, controlPts] = convert3DNURBS(solid);

    [element, index, noElems, ...
     elRangeXi, elRangeEta, elRangeZeta] = generateIGA3DMesh(Xi, Eta, Zeta, p, q, r, l, m, n);

    % build visualization B8 mesh
    vtfFileName = ['../graphics/GLview/rectangularPrismAnalytic' num2str(analyticFunction)];

    noUniqueXiKnots = length(unique(Xi));
    noUniqueEtaKnots = length(unique(Eta));
    noUniqueZetaKnots = length(unique(Zeta));

    extraXiPts = floor(20/(noUniqueXiKnots-1)); % .. per element
    extraEtaPts = floor(20/(noUniqueEtaKnots-1)); % .. per element
    extraZetaPts = floor(20/(noUniqueZetaKnots-1)); % .. per element
    % 
    [nodes, noNodes, visElements, cornerNode, ...
              noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, solid);
    % 
    displacement   = zeros(noNodes,3);
    stress   = zeros(noNodes,6);

    for i = 1:noNodes
        v = nodes(i,1:3);
        x = v(1);
        y = v(2);
        z = v(3);

        displacement(i,:) = u_analytic(x,y,z)';

        stress(i,:) = stress_u(x,y,z)';
    end

    dataAnal.nodes = nodes;
    dataAnal.visElements = visElements;
    dataAnal.displacement = displacement;
    dataAnal.stress = stress;
    options = {'name',vtfFileName,...
         'plotDisplacementVectors',1,...
         'plotVonMisesStress',1};

    makeVTFfile_new(dataAnal, options);
end

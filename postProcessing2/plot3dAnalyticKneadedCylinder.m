close all
clear all
for analyticFunction = 1:3
    analyticFunction
    switch analyticFunction
        case 1
            c_1 = 1;
            c_2 = 1;
            c_3 = 1;
        case 2
            c_1 = 7;
            c_2 = 5;
            c_3 = 1;
        case 3
            c_1 = 4;
            c_2 = 2;
            c_3 = 1;        
    end
    R = 1.5; % Mid surface radius
    t = 1;  % Thickness
    R_i = R-t/2;    % Inner surface radius
    R_o = R+t/2;    % Outer surface radius
    L = 5;  % Height
    p_i = 1;
    p_o = 0;
    E = 13;
    nu = 0.3;
    lambda = nu*E/(1+nu)/(1-2*nu);
    mu = E/(2*(1+nu));

    getAnalyticFunctions_kneadedCylinder

    solid = getSolidCylinderData(R_i, R_o, L);

    [Xi, Eta, Zeta, ...
     p, q, r, ...
     n, m, l, ...
     weights, controlPts] = convert3DNURBS(solid);

    [element, index, noElems, ...
     elRangeXi, elRangeEta, elRangeZeta] = generateIGA3DMesh(Xi, Eta, Zeta, p, q, r, l, m, n);

    % build visualization B8 mesh
    vtfFileName = ['../graphics/GLview/kneadedCylinderAnalytic' num2str(analyticFunction)];

    noUniqueXiKnots = length(unique(Xi));
    noUniqueEtaKnots = length(unique(Eta));
    noUniqueZetaKnots = length(unique(Zeta));

    extraXiPts = floor(100/(noUniqueXiKnots-1)); % .. per element
    extraEtaPts = floor(100/(noUniqueEtaKnots-1)); % .. per element
    extraZetaPts = floor(0/(noUniqueZetaKnots-1)); % .. per element
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
        radius = sqrt(x^2+y^2);
        theta = atan2(y,x);

        displacement(i,:) = u_analytic(radius,theta,z)';

        stress(i,:) = stress_u(radius,theta,z)';
    end

    dataAnal.nodes = nodes;
    dataAnal.visElements = visElements;
    dataAnal.displacement = displacement;
    dataAnal.stress = stress;
    options = {'name',vtfFileName,...
         'plotDisplacementVectors',1,...
         'plotPolarStress_thetaz',1,...
         'plotMagnitudeDisplacement',1,...
         'plotPolarRadialDisplacement',1};

    makeVTFfile_new(dataAnal, options);
end

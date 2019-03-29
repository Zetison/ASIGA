close all
clear all

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

getAnalyticFunctions_solidCircularCylinder

solid = getSolidCylinderData(R_i, R_o, L);

[Xi, Eta, Zeta, ...
 p, q, r, ...
 n, m, l, ...
 weights, controlPts] = convert3DNURBS(solid);

[element, index, noElems, ...
 elRangeXi, elRangeEta, elRangeZeta] = generateIGA3DMesh(Xi, Eta, Zeta, p, q, r, l, m, n);

% build visualization B8 mesh
vtfFileName = '../graphics/GLview/solidCylinderAnalytic';
            
noUniqueXiKnots = length(unique(Xi));
noUniqueEtaKnots = length(unique(Eta));
noUniqueZetaKnots = length(unique(Zeta));

extraXiPts = floor(500/(noUniqueXiKnots-1)); % .. per element
extraEtaPts = floor(0/(noUniqueEtaKnots-1)); % .. per element
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
    radius = sqrt(x^2+y^2);
    theta = atan2(y,x);
    u_r_anal = u_analytic(radius,theta,z);
    
    displacement(i,1) = u_r_anal*cos(theta);
    displacement(i,2) = u_r_anal*sin(theta);
    displacement(i,3) = 0;
    
    stress(i,:) = stress_u(radius,theta,z)';
end

dataAnal.nodes = nodes;
dataAnal.visElements = visElements;
dataAnal.displacement = displacement;
dataAnal.stress = stress;
options = {'name',vtfFileName,...
     'plotDisplacementVectors',1,...
     'plotPolarRadialDisplacement',1};


makeVTFfile_new(dataAnal, options);


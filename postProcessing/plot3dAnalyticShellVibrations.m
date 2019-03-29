% getSolidCylinderData
% 
% convert3DNURBS

% build visualization B8 mesh
%             
noUniqueXiKnots = length(unique(Xi));
noUniqueEtaKnots = length(unique(Eta));
noUniqueZetaKnots = length(unique(Zeta));
% 
extraXiPts = floor(50/(noUniqueXiKnots-1)); % .. per element
extraEtaPts = floor(50/(noUniqueEtaKnots-1)); % .. per element
extraZetaPts = floor(0/(noUniqueZetaKnots-1)); % .. per element

[nodes, noNodes, visElements, cornerNode, ...
 noXiKnots, noEtaKnots, noZetaKnots] = buildVisualization3dMesh_new(Xi, Eta, Zeta, extraXiPts, extraEtaPts, extraZetaPts, noElems, solid);


[~, analyticEigenValues] = getDataFromFile2('plotData/sphericalShell/eigenvalues.dat');
for vibMode = 0:10
    vtfFileName = ['../graphics/GLview/vibrationsInVacuoSphericalShellAnalytic' num2str(vibMode)];
    omega = analyticEigenValues(vibMode+1);
    displacement   = zeros(noNodes,3);

    n = vibMode
    lambda = nu*E/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));
    c_1 = sqrt((lambda+2*mu)/rho_s);
    c_2 = sqrt(mu/rho_s);

    alpha = omega/c_1;
    beta = omega/c_2;

    for j = 1:noNodes
        v = nodes(j,1:3);
        x = v(1);
        y = v(2);
        z = v(3);
        radius = sqrt(x^2+y^2+z^2);
        phi = atan2(y,x);
        theta = acos(z/radius);
        
        u_r = 1/radius*(U_fun(1,1,n,alpha*radius)+U_fun(3,1,n,beta*radius))*legendre3(n,cos(theta)) ...
             +1/radius*(U_fun(1,2,n,alpha*radius)+U_fun(3,2,n,beta*radius))*legendre3(n,cos(theta));
        if abs(theta) < newEpsilon || abs(theta - pi) < newEpsilon
            u_theta = 0;
        else
            u_theta = 1/radius*(V_fun(1,1,n,alpha*radius) + V_fun(3,1,n,beta*radius) ...
                                 + V_fun(1,2,n,alpha*radius) + V_fun(3,2,n,beta*radius)) ...
                           *n*(cot(theta)*legendre3(n,cos(theta))-1/sin(theta)*legendre3(n-1,cos(theta)));
        end
        u_phi = 0;


        displacement(j,1) = u_r*cos(phi)*sin(theta) + u_phi*(-sin(phi)) + u_theta*cos(phi)*cos(theta);
        displacement(j,2) = u_r*sin(phi)*sin(theta) + u_phi*cos(phi)    + u_theta*sin(phi)*cos(theta);
        displacement(j,3) = u_r*cos(theta)          + u_phi*0           + u_theta*(-sin(theta));

        if isnan(displacement(j,1)) || isnan(displacement(j,2)) || isnan(displacement(j,3))
            keyboard;
        end
    %     u_r = 2; %sin(7*pi*z/H);
    %     u_theta = 0; %sin(7*pi*z/H);
    %     u_z = 0; %sin(7*pi*z/H);
    %     
    %     displacement(i,1) = u_r*cos(theta)-u_theta*sin(theta);
    %     displacement(i,2) = u_r*sin(theta)+u_theta*cos(theta);
    %     displacement(i,3) = u_z;
    end
    dataAnal.omega = omega;
    dataAnal.nodes = nodes;
    dataAnal.visElements = visElements;
    dataAnal.displacement = displacement;
    options = {'name',vtfFileName,...
         'plotDisplacementVectors',1,...
         'plotSphericalRadialDisplacement',1,...
         'plotTimeOscillation',0};


    makeVTFfile_new(dataAnal, options);
end


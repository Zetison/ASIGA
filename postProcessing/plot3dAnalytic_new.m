% getSolidCylinderData
% 
% convert3DNURBS

% build visualization B8 mesh
vtfFileName = '../graphics/GLview/kneadedCylinderAnalytic';
%             
% noUniqueXiKnots = length(unique(Xi));
% noUniqueEtaKnots = length(unique(Eta));
% noUniqueZetaKnots = length(unique(Zeta));
% 
% extraXiPts = floor(100/(noUniqueXiKnots-1)); % .. per element
% extraEtaPts = floor(100/(noUniqueEtaKnots-1)); % .. per element
% extraZetaPts = floor(7/(noUniqueZetaKnots-1)); % .. per element
% 
% buildVisualization3dMesh_new
% 
% displacement   = zeros(noNodes,3);

for i = 1:noNodes
    v = nodes(i,1:3);
    x = v(1);
    y = v(2);
    z = v(3);
    radius = sqrt(x^2+y^2);
    theta = atan2(y,x);
    u_anal = u_analytic(radius,theta,z);
    
    displacement(i,1) = u_anal(1);
    displacement(i,2) = u_anal(2);
    displacement(i,3) = u_anal(3);
    
%     u_r = 2; %sin(7*pi*z/H);
%     u_theta = 0; %sin(7*pi*z/H);
%     u_z = 0; %sin(7*pi*z/H);
%     
%     displacement(i,1) = u_r*cos(theta)-u_theta*sin(theta);
%     displacement(i,2) = u_r*sin(theta)+u_theta*cos(theta);
%     displacement(i,3) = u_z;
end

dataAnal.nodes = nodes;
dataAnal.visElements = visElements;
dataAnal.displacement = displacement;
options = {'name',vtfFileName,...
     'plotDisplacementVectors',1,...
     'plotXdisplacement',1,...
     'plotYdisplacement',1,...
     'plotZdisplacement',1,...
     'plotPolarRadialDisplacement',1,...
     'plotStressXX',0,...
     'plotStressYY',0,...
     'plotStressZZ',0,...
     'plotStressYZ',0,...
     'plotStressXZ',0,...
     'plotStressXY',0,...
     'plotVonMisesStress',0};


makeVTFfile_new(dataAnal, options);


function makeVTFfile(data, newOptions)

options = struct('name','untitled',...
                 'plotDisplacementVectors',0,...
                 'plotScalarField',0,...
                 'plotPolarRadialDisplacement',0,...
                 'plotSphericalRadialDisplacement',0,...
                 'plotStressXX',0,...
                 'plotStressYY',0,...
                 'plotStressZZ',0,...
                 'plotStressYZ',0,...
                 'plotStressXZ',0,...
                 'plotStressXY',0,...
                 'plotPolarStress_rr',0,...
                 'plotPolarStress_thetatheta',0,...
                 'plotPolarStress_zz',0,...
                 'plotPolarStress_thetaz',0,...
                 'plotPolarStress_rz',0,...
                 'plotPolarStress_rtheta',0,...
                 'plotSphericalStress_rr',0,...
                 'plotSphericalStress_thetatheta',0,...
                 'plotSphericalStress_phiphi',0,...
                 'plotSphericalStress_thetaphi',0,...
                 'plotSphericalStress_rphi',0,...
                 'plotSphericalStress_rtheta',0,...
                 'plotVonMisesStress',0,...
                 'plotOrgJacobi',0,...
                 'plotTimeOscillation', 0,...
                 'plotError', 0,...
                 'plotErrorGradient', 0,...
                 'plotErrorTot', 0,...
                 'plotGradientVectors', 0, ...
                 'plotP_inc', 0, ...
                 'plotTotField',0,...
                 'addDynamicScalars', 0);

%# read the acceptable names
optionNames = fieldnames(options);

%# count arguments
nArgs = length(newOptions);
if round(nArgs/2) ~= nArgs/2
   error('Must have propertyName/propertyValue pairs')
end

for pair = reshape(newOptions,2,[]) %# pair is {propName;propValue}
   inpName = pair{1}; %# make case insensitive

   if any(strcmp(inpName,optionNames))
      %# overwrite options. If you want you can test for the right class here
      %# Also, if you find out that there is an option you keep getting wrong,
      %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
      options.(inpName) = pair{2};
   else
      error('%s is not a recognized parameter name',inpName)
   end
end




vtfFileName = options.name;
nodes = data.nodes;
visElements = data.visElements;
if options.plotTimeOscillation
    noSteps = 30;
    omega = data.omega;
else
    noSteps = 1;
    omega = 1;
end

%--------------------------------------------------------------------------
% Write the header for a structured grid:
%--------------------------------------------------------------------------
fid = fopen([vtfFileName,'.vtf'],'wt+','b'); % note the 'b': not doing BigEndian does not work with MATLAB!
fprintf(fid,'*VTF-1.00\n');

noVisElems = size(visElements,1);
noNodes = size(nodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add coordinates
fprintf(fid,'\n*NODES 1\n');
fprintf(fid,'%%WITH_ID \n');
for i = 1:noNodes
    fprintf(fid,' %d %21.15f %21.15f %21.15f\n', [i nodes(i,1) nodes(i,2) nodes(i,3)]);
end

fprintf(fid,'\n*ELEMENTS 1\n');
fprintf(fid,'%%NAME "Hex elements" \n');
fprintf(fid,'%%DESCRIPTION "The hexahedron elements" \n');
fprintf(fid,'%%NODES #1 \n');
fprintf(fid,'%%WITH_ID \n');
fprintf(fid,'%%HEXAHEDRONS \n');
for eVis = 1:noVisElems
    fprintf(fid,' %d %d %d %d %d %d %d %d %d\n', [eVis visElements(eVis,:)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add results
counter = 1;
%% Displacement (vectors)
if options.plotDisplacementVectors
    counter = printField(data.displacement, nodes, noSteps, omega, counter, fid);
end

%% Gradient (vectors)
if options.plotGradientVectors
    counter = printField(data.gradient, nodes, noSteps, omega, counter, fid);
end

%% Scalar field
if options.plotScalarField
    counter = printField(data.scalarField, nodes, noSteps, omega, counter, fid);
end

%% Total scalar field
if options.plotTotField
    counter = printField(data.totField, nodes, noSteps, omega, counter, fid);
end

%% P_inc
if options.plotP_inc
    counter = printField(data.P_inc, nodes, noSteps, omega, counter, fid);
end

%% Error in scalar field
if options.plotError
    counter = printField(data.Error, nodes, 1, omega, counter, fid);
end

%% Error in gradient field
if options.plotErrorGradient
    counter = printField(data.ErrorGradient, nodes, 1, omega, counter, fid);
end

%% Total Error
if options.plotErrorTot
    counter = printField(data.errorTot, nodes, 1, omega, counter, fid);
end

%% Displacement polar radial direction (scalars)
if options.plotPolarRadialDisplacement
    x = nodes(:,1);
    y = nodes(:,2);
    theta = atan2(y,x);
    displacement = data.displacement;
    u_x = displacement(:,1);
    u_y = displacement(:,2);
    u_r = u_y.*sin(theta) + u_x.*cos(theta);
    counter = printField(u_r, nodes, noSteps, omega, counter, fid);
end

%% Displacement spherical radial direction (scalars)
if options.plotPolarRadialDisplacement
    x = nodes(:,1);
    y = nodes(:,2);
    z = nodes(:,3);
    r = sqrt(x^2+y^2+z^2);
    phi = atan2(y,x);
    theta = acos(z./r);
    displacement = data.displacement;
    u_x = displacement(:,1);
    u_y = displacement(:,2);
    u_z = displacement(:,3);
    u_r = u_x*sin(theta)*cos(phi) + u_y*sin(theta)*sin(phi) + u_z*cos(theta);
%     u_theta = u_x*cos(theta)*cos(phi) + u_y*cos(theta)*sin(phi) - u_z*sin(theta);
%     u_phi = -u_x*sin(phi) + u_y*cos(phi);
    counter = printField(u_r, nodes, noSteps, omega, counter, fid);
end

%% Stress (scalars)
toggleStressVector = [options.plotStressXX options.plotStressYY options.plotStressZZ options.plotStressYZ options.plotStressXZ options.plotStressXY];
for i = 1:6
    if toggleStressVector(i)
        counter = printField(data.stress(:,i), nodes, noSteps, omega, counter, fid);
    end
end

%% Polar stress (scalars)
togglePolarStressVector = [options.plotPolarStress_rr options.plotPolarStress_thetatheta options.plotPolarStress_zz options.plotPolarStress_thetaz options.plotPolarStress_rz options.plotPolarStress_rtheta];
for i = 1:6
    if togglePolarStressVector(i)
        x = nodes(:,1);
        y = nodes(:,2);
        theta = atan2(y,x);
        stress = data.stress;
        polarStress = [(cos(theta))^2*stress(:,1) + (sin(theta))^2*stress(:,2) + sin(2*theta)*stress(:,6), ...
                       (sin(theta))^2*stress(:,1) + (cos(theta))^2*stress(:,2) - sin(2*theta)*stress(:,6), ...
                       stress(:,3), ...
                       stress(:,4)*cos(theta) - stress(:,5)*sin(theta), ...
                       stress(:,4)*sin(theta) + stress(:,5)*cos(theta), ...
                       -1/2*sin(2*theta)*stress(:,1) + 1/2*sin(2*theta)*stress(:,2) - cos(2*theta)*stress(:,6)];
        
        counter = printField(polarStress(:,i), nodes, noSteps, omega, counter, fid);
    end
end

%% Spherical stress (scalars) 
toggleSphericalStressVector = [options.plotSphericalStress_rr options.plotSphericalStress_thetatheta options.plotSphericalStress_phiphi ...
                                    options.plotSphericalStress_thetaphi options.plotSphericalStress_rphi options.plotSphericalStress_rtheta];
for i = 1:6
    if toggleSphericalStressVector(i)
        x = nodes(:,1);
        y = nodes(:,2);
        z = nodes(:,3);
        r = sqrt(x.^2+y.^2+z.^2);
        phi_arr = atan2(y,x);
        theta_arr = acos(z./r);
        stress = data.stress;
        sphericalStress = zeros(size(stress));
        for j = 1:size(theta_arr,1)
            theta = theta_arr(j);
            phi = phi_arr(j);
            D = [sin(theta)^2*cos(phi)^2, sin(theta)^2*sin(phi)^2, cos(theta)^2,  sin(2*theta)*sin(phi),  sin(2*theta)*cos(phi), sin(theta)^2*sin(2*phi);
                 cos(theta)^2*cos(phi)^2, cos(theta)^2*sin(phi)^2, sin(theta)^2, -sin(2*theta)*sin(phi), -sin(2*theta)*cos(phi), cos(theta)^2*sin(2*phi);
                 sin(phi)^2,              cos(phi)^2,              0,             0,                      0,                     -sin(2*phi);
                 -0.5*cos(theta)*sin(2*phi), 0.5*cos(theta)*sin(2*phi), 0,        -sin(theta)*cos(phi),   sin(theta)*sin(phi),   cos(theta)*cos(2*phi);
                 -0.5*sin(theta)*sin(2*phi), 0.5*sin(theta)*sin(2*phi), 0,        cos(theta)*cos(phi),    -cos(theta)*sin(phi),   sin(theta)*cos(2*phi);
                  0.5*sin(2*theta)*cos(phi)^2, 0.5*sin(2*theta)*sin(phi)^2, -0.5*sin(2*theta),        cos(2*theta)*sin(phi),    cos(2*theta)*cos(phi),   0.5*sin(2*theta)*sin(2*phi)];

            sphericalStress(j,:) = (D*stress(j,:).').';
        end
        
        counter = printField(sphericalStress(:,i), nodes, noSteps, omega, counter, fid);
    end
end

%% von Mises Stress (scalars)
if options.plotVonMisesStress
    stress = data.stress;
    sigma_v = sqrt((stress(:,1)-stress(:,2)).^2 + (stress(:,2)-stress(:,3)).^2 + (stress(:,1)-stress(:,3)).^2 ...
                                + 6*(stress(:,4).^2 + stress(:,5).^2 + stress(:,6).^2)/2);
    counter = printField(sigma_v, nodes, noSteps, omega, counter, fid);
end

%% Jacobian (scalars)
if options.plotOrgJacobi
    printField(data.orgJacobi, nodes, noSteps, omega, counter, fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine results

%% Geometry
fprintf(fid,'\n*GLVIEWGEOMETRY 1\n');
fprintf(fid,'%%ELEMENTS\n');
fprintf(fid,'1\n');

glviewScalarCounter = 1;
glviewVectorCounter = 1;
counter = 1;

%% Displacement (vector)
if options.plotDisplacementVectors
    [counter, glviewVectorCounter] = printGeometry(counter, glviewVectorCounter, noSteps, 3, 'Displacement', fid);
end

%% Gradient (vector)
if options.plotGradientVectors
    counter = printGeometry(counter, glviewVectorCounter, noSteps, 3, 'Gradient', fid);
end

%% Scalar field (scalars)
if options.plotScalarField
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Scalar field', fid);
end

%% Total scalar field (scalars)
if options.plotTotField
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Total scalar field', fid);
end

%% P_inc (scalars)
if options.plotP_inc
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'P_inc', fid);
end

%% Error in scalar field (scalars)
if options.plotError
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Error', fid);
end

%% Error in gradient field (scalars)
if options.plotErrorGradient
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Error gradient', fid);
end

%% Total error (scalars)
if options.plotErrorTot
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Total error', fid);
end


%% Displacement polar radial direction (scalars)
if options.plotPolarRadialDisplacement
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Displacement polar radial direction', fid);
end

%% Displacement spherical radial direction (scalars)
if options.plotSphericalRadialDisplacement
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Displacement spherical radial direction', fid);
end

%% Stress (scalars)
stressTitles = {'Stress xx','Stress yy','Stress zz','Stress yz','Stress xz','Stress xy'};
for i = 1:6
    if toggleStressVector(i)
        [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, stressTitles{i}, fid);
    end
end

%% Polar stress (scalars)
stressTitles = {'Stress rr','Stress thetatheta','Stress zz','Stress thetaz','Stress rz','Stress rtheta'};
for i = 1:6
    if togglePolarStressVector(i)
        [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, stressTitles{i}, fid);
    end
end

%% Spherical stress (scalars)
stressTitles = {'Stress rr','Stress thetatheta','Stress phiphi','Stress thetaphi','Stress rphi','Stress rtheta'};
for i = 1:6
    if toggleSphericalStressVector(i)
        [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, stressTitles{i}, fid);
    end
end

%% von Mises Stress (scalars)
if options.plotVonMisesStress
    [counter, glviewScalarCounter] = printGeometry(counter, glviewScalarCounter, noSteps, 1, 'von Mises Stress', fid);
end

%% Jacobian (scalars)
if options.plotOrgJacobi
    printGeometry(counter, glviewScalarCounter, noSteps, 1, 'Jacobian', fid);
end

fclose(fid);

function resultCounter = printField(field, nodes, noSteps, omega, resultCounter, fid)
d = size(field,2);
for s = 1:noSteps
    t = (s-1)/noSteps*2*pi/omega;

    fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
    fprintf(fid,'%%WITH_ID \n');
    fprintf(fid,'%%DIMENSION %d\n', d);
    fprintf(fid,'%%PER_NODE #1\n');
    for i = 1:size(nodes,1)
        fprintf(fid,[' %d' repmat(' %21.15f',1,d) '\n'], [i real(field(i,:)*exp(-1i*omega*t))]);
    end
    resultCounter = resultCounter + 1;
end

function [counter, counter2] = printGeometry(counter, counter2, noSteps, d, name, fid)
if d == 1
    fprintf(fid,'\n*GLVIEWSCALAR %d\n', counter2);
else
    fprintf(fid,'\n*GLVIEWVECTOR %d\n', counter2);
end
counter2 = counter2 + 1;

fprintf(fid,'%%NAME "%s"\n', name);
for s = 1:noSteps
    fprintf(fid,'%%STEP %d\n', s);
    fprintf(fid,'%d\n', counter);
    counter = counter + 1;
end


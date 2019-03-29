function makeVisualizationFile(data, newOptions)

options = struct('name','untitled',...
                 'plotDisplacementVectors',0,...
                 'plotScalarField',0,...
                 'plotXdisplacement',0,...
                 'plotYdisplacement',0,...
                 'plotZdisplacement',0,...
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
                 'plotVonMisesStress',0,...
                 'plotOrgJacobi',0,...
                 'plotMagnitudeDisplacement',0,...
                 'plotTimeOscillation', 0,...
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
noVisElems = size(visElements,1);
noNodes = size(nodes,1);
switch options.filetype
    case 'vtf'
        fid = fopen([vtfFileName,'.vtf'],'wt+','b');
        fprintf(fid,'*VTF-1.00\n');
    case 'vtk'
        fid = fopen([vtfFileName,'.vtk'],'wt+','b');
        fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
        fprintf(fid,'\t<UnstructuredGrid>\n');
        fprintf(fid,'\t\t<Piece NumberOfPoints="%d" NumberOfCells="%d">',noNodes,noVisElems);
end
                                
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add coordinates
switch options.filetype
    case 'vtf'
        fprintf(fid,'\n*NODES 1\n');
        fprintf(fid,'%%WITH_ID \n');
        for i = 1:noNodes
            fprintf(fid,' %d %g %g %g\n', [nodes(i,4) nodes(i,1) nodes(i,2) nodes(i,3)]);
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
    case 'vtk'
        fprintf(fid,'\t\t\t<Points>\n');
        fprintf(fid,'\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">');
        for i = 1:noNodes
            fprintf(fid,'\t\t\t\t\t\t%g %g %g\n', [nodes(i,1) nodes(i,2) nodes(i,3)]);
        end
        fprintf(fid,'\t\t\t\t</DataArray>\n');
        fprintf(fid,'\t\t\t</Points>\n');
        
        
        fprintf(fid,'\t\t\t<Cells>\n');
        fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n');
        for eVis = 1:noVisElems
            fprintf(fid,'\t\t\t\t\t\t%d %d %d %d %d %d %d %d\n', visElements(eVis,:));
        end
        fprintf(fid,'\t\t\t\t</DataArray>\n');
        
        fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n');
        fprintf(fid,'\t\t\t\t\t\t8 16 24 32 40 48 56 64\n');
        fprintf(fid,'\t\t\t\t</DataArray>\n');
        
        fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="types" format="ascii">\n');
        fprintf(fid,'\t\t\t\t\t\t11 11 11 11 11 11 11 11\n');
        fprintf(fid,'\t\t\t\t</DataArray>\n');
        fprintf(fid,'\t\t\t</Cells>\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add results
switch options.filetype
    case 'vtf'
        resultCounter = 1;
    case 'vtk'
        fprintf(fid,'\t\t\t<PointData Scalars="scalars">\n');
end
       
%% Displacement (vectors)
if options.plotDisplacementVectors    
    switch options.filetype
        case 'vtf'
            for s = 1:noSteps
                t = (s-1)/noSteps*2*pi/omega;
                displacement = data.displacement;

                fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
                fprintf(fid,'%%WITH_ID \n');
                fprintf(fid,'%%DIMENSION 3 \n');
                fprintf(fid,'%%PER_NODE #1 \n');
                for i = 1:noNodes
                    fprintf(fid,' %d %g %g %g\n', [nodes(i,4) displacement(i,:)*cos(omega*t)]);
                end

                resultCounter = resultCounter + 1;
            end
        case 'vtk'
            fprintf(fid,'\t\t\t<Points>\n');
            fprintf(fid,'\t\t\t\t<DataArray type="Float32" NumberOfComponents="3" format="ascii">');
            for i = 1:noNodes
                fprintf(fid,'\t\t\t\t\t\t%g %g %g\n', [nodes(i,1) nodes(i,2) nodes(i,3)]);
            end
            fprintf(fid,'\t\t\t\t</DataArray>\n');
            fprintf(fid,'\t\t\t</Points>\n');
    end
            
end

%% Scalar field
if options.plotScalarField
    if options.addDynamicScalars
        noStepsUpdated = noSteps;
    else
        noStepsUpdated = 1;
    end
    for s = 1:noStepsUpdated
        t = (s-1)/noStepsUpdated*2*pi/omega;
        scalarField = data.scalarField;
        fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
        fprintf(fid,'%%WITH_ID \n');
        fprintf(fid,'%%DIMENSION 1 \n');
        fprintf(fid,'%%PER_NODE #1 \n');
        for i = 1:noNodes
            fprintf(fid,' %d %g\n', [nodes(i,4) scalarField(i)*cos(omega*t)]);
        end
    resultCounter = resultCounter + 1;
    end
end

%% Displacement (scalars)
toggleDispVector = [options.plotXdisplacement options.plotYdisplacement options.plotZdisplacement];
for d = 1:3
    if toggleDispVector(d)
        if options.addDynamicScalars
            noStepsUpdated = noSteps;
        else
            noStepsUpdated = 1;
        end
        for s = 1:noStepsUpdated
            t = (s-1)/noStepsUpdated*2*pi/omega;
            displacement = data.displacement;
            fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
            fprintf(fid,'%%WITH_ID \n');
            fprintf(fid,'%%DIMENSION 1 \n');
            fprintf(fid,'%%PER_NODE #1 \n');
            for i = 1:noNodes
                fprintf(fid,' %d %g\n', [nodes(i,4) displacement(i,d)*cos(omega*t)]);
            end
        resultCounter = resultCounter + 1;
        end
    end
end

%% Displacement magnitude (scalars)
if options.plotMagnitudeDisplacement
    if options.addDynamicScalars
        noStepsUpdated = noSteps;
    else
        noStepsUpdated = 1;
    end
    for s = 1:noStepsUpdated
        t = (s-1)/noStepsUpdated*2*pi/omega;
        displacement = data.displacement;
        fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
        fprintf(fid,'%%WITH_ID \n');
        fprintf(fid,'%%DIMENSION 1 \n');
        fprintf(fid,'%%PER_NODE #1 \n');
        for i = 1:noNodes
            u = displacement(i,:);
            fprintf(fid,' %d %g\n', [nodes(i,4) norm(u,2)*cos(omega*t)]);
        end
        resultCounter = resultCounter + 1;
    end
end

%% Displacement polar radial direction (scalars)
if options.plotPolarRadialDisplacement
    if options.addDynamicScalars
        noStepsUpdated = noSteps;
    else
        noStepsUpdated = 1;
    end
    for s = 1:noStepsUpdated
        t = (s-1)/noStepsUpdated*2*pi/omega;
        displacement = data.displacement;
        fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
        fprintf(fid,'%%WITH_ID \n');
        fprintf(fid,'%%DIMENSION 1 \n');
        fprintf(fid,'%%PER_NODE #1 \n');
        for i = 1:noNodes
            x = nodes(i,1);
            y = nodes(i,2);
            u_x = displacement(i,1);
            u_y = displacement(i,2);
            theta = atan2(y,x);
            u_r = u_y*sin(theta) + u_x*cos(theta);
            fprintf(fid,' %d %g\n', [nodes(i,4) u_r*cos(omega*t)]);
        end
        resultCounter = resultCounter + 1;
    end
end

%% Displacement spherical radial direction (scalars)
if options.plotSphericalRadialDisplacement
    if options.addDynamicScalars
        noStepsUpdated = noSteps;
    else
        noStepsUpdated = 1;
    end
    for s = 1:noStepsUpdated
        t = (s-1)/noStepsUpdated*2*pi/omega;
        displacement = data.displacement;
        fprintf(fid,['\n*RESULTS ' num2str(resultCounter) '\n']);
        fprintf(fid,'%%WITH_ID \n');
        fprintf(fid,'%%DIMENSION 1 \n');
        fprintf(fid,'%%PER_NODE #1 \n');
        for i = 1:noNodes
            x = nodes(i,1);
            y = nodes(i,2);
            z = nodes(i,3);
            r = sqrt(x^2+y^2+z^2);
            phi = atan2(y,x);
            theta = acos(z/r);
            
            u_x = displacement(i,1);
            u_y = displacement(i,2);
            u_z = displacement(i,3);
            
            u_r = u_x*sin(theta)*cos(phi) + u_y*sin(theta)*sin(phi) + u_z*cos(theta);
%             u_theta = u_x*cos(theta)*cos(phi) + u_y*cos(theta)*sin(phi) - u_z*sin(theta);
%             u_phi = -u_x*sin(phi) + u_y*cos(phi);
            
            fprintf(fid,' %d %g\n', [nodes(i,4) u_r*cos(omega*t)]);
        end
        resultCounter = resultCounter + 1;
    end
end

%% Stress (scalars)
toggleStressVector = [options.plotStressXX options.plotStressYY options.plotStressZZ options.plotStressYZ options.plotStressXZ options.plotStressXY];
for d = 1:6
    if toggleStressVector(d)
        stress = data.stress;
        if options.addDynamicScalars
            noStepsUpdated = noSteps;
        else
            noStepsUpdated = 1;
        end
        for s = 1:noStepsUpdated
            t = (s-1)/noStepsUpdated*2*pi/omega;
            fprintf(fid,['\n*RESULTS ' num2str(resultCounter) ' \n']);
            fprintf(fid,'%%WITH_ID \n');
            fprintf(fid,'%%DIMENSION 1 \n');
            fprintf(fid,'%%PER_NODE #1 \n');
            for i = 1:noNodes
                fprintf(fid,' %d %g\n', [nodes(i,4) stress(i,d)*cos(omega*t)]);
            end
            resultCounter = resultCounter + 1;
        end
    end
end

%% Polar stress (scalars)
togglePolarStressVector = [options.plotPolarStress_rr options.plotPolarStress_thetatheta options.plotPolarStress_zz options.plotPolarStress_thetaz options.plotPolarStress_rz options.plotPolarStress_rtheta];
for d = 1:6
    if togglePolarStressVector(d)
        stress = data.stress;
        polarStress = zeros(size(stress));
        for i = 1:size(data.stress,1)
            x = nodes(i,1);
            y = nodes(i,2);
            theta = atan2(y,x);
            
            polarStress(i,:) = [(cos(theta))^2*stress(i,1) + (sin(theta))^2*stress(i,2) + sin(2*theta)*stress(i,6);
                           (sin(theta))^2*stress(i,1) + (cos(theta))^2*stress(i,2) - sin(2*theta)*stress(i,6);
                           stress(i,3);
                           stress(i,4)*cos(theta) - stress(i,5)*sin(theta);
                           stress(i,4)*sin(theta) + stress(i,5)*cos(theta);
                           -1/2*sin(2*theta)*stress(i,1) + 1/2*sin(2*theta)*stress(i,2) - cos(2*theta)*stress(i,6)]';
        end
        if options.addDynamicScalars
            noStepsUpdated = noSteps;
        else
            noStepsUpdated = 1;
        end
        for s = 1:noStepsUpdated
            t = (s-1)/noStepsUpdated*2*pi/omega;
            fprintf(fid,['\n*RESULTS ' num2str(resultCounter) ' \n']);
            fprintf(fid,'%%WITH_ID \n');
            fprintf(fid,'%%DIMENSION 1 \n');
            fprintf(fid,'%%PER_NODE #1 \n');
            for i = 1:noNodes
                fprintf(fid,' %d %g\n', [nodes(i,4) polarStress(i,d)*cos(omega*t)]);
            end
            resultCounter = resultCounter + 1;
        end
    end
end

%% von Mises Stress (scalars)
if options.plotVonMisesStress
	stress = data.stress;
    if options.addDynamicScalars
        noStepsUpdated = noSteps;
    else
        noStepsUpdated = 1;
    end
    for s = 1:noStepsUpdated
        t = (s-1)/noStepsUpdated*2*pi/omega;
        fprintf(fid,['\n*RESULTS ' num2str(resultCounter) ' \n']);
        fprintf(fid,'%%WITH_ID \n');
        fprintf(fid,'%%DIMENSION 1 \n');
        fprintf(fid,'%%PER_NODE #1 \n');
        for i = 1:noNodes
            sigma_v = sqrt((stress(i,1)-stress(i,2))^2 + (stress(i,2)-stress(i,3))^2 + (stress(i,1)-stress(i,3))^2 ...
                                        + 6*(stress(i,4)^2 + stress(i,5)^2 + stress(i,6)^2)/2);
            fprintf(fid,' %d %g\n', [nodes(i,4) sigma_v*abs(cos(omega*t))]);
        end 
        resultCounter = resultCounter + 1;
    end
end

%% Jacobian (scalars)
if options.plotOrgJacobi
	orgJacobi = data.orgJacobi;
    
    if options.addDynamicScalars
        noStepsUpdated = noSteps;
    else
        noStepsUpdated = 1;
    end
    for s = 1:noStepsUpdated
        fprintf(fid,['\n*RESULTS ' num2str(resultCounter) ' \n']);
        fprintf(fid,'%%WITH_ID \n');
        fprintf(fid,'%%DIMENSION 1 \n');
        fprintf(fid,'%%PER_NODE #1 \n');
        for i = 1:noNodes
            fprintf(fid,' %d %g\n', [nodes(i,4) orgJacobi(i)]);
        end 
        resultCounter = resultCounter + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine results

switch options.filetype
    case 'vtf'
        %% Geometry
        fprintf(fid,'\n*GLVIEWGEOMETRY 1\n');
        fprintf(fid,'%%ELEMENTS\n');
        fprintf(fid,'1\n');        

        %% Displacement (vector)
        resultCounter = 1;
        if options.plotDisplacementVectors
            fprintf(fid,'\n*GLVIEWVECTOR 1\n');
            fprintf(fid,'%%NAME "Displacement"\n');
            for s = 1:noSteps
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            resultCounter = resultCounter + noSteps;
        end

        %% Scalar field (scalars)
        glviewScalarCounter = 1;
        if options.plotScalarField
            fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
            fprintf(fid,'%%NAME "Scalar field"\n');
            if options.addDynamicScalars
                noStepsUpdated = noSteps;
            else
                noStepsUpdated = 1;
            end
            for s = 1:noStepsUpdated
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            if options.addDynamicScalars
                resultCounter = resultCounter + noStepsUpdated;
            else
                resultCounter = resultCounter + 1;
            end
            glviewScalarCounter = glviewScalarCounter + 1;
        end

        %% Displacement (scalars)
        dispTitles = {'x displacement','y displacement','z displacement'};
        for d = 1:3
            if toggleDispVector(d)
                fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
                fprintf(fid,['%%NAME "' dispTitles{d} '"\n']);
                if options.addDynamicScalars
                    noStepsUpdated = noSteps;
                else
                    noStepsUpdated = 1;
                end
                for s = 1:noStepsUpdated
                    fprintf(fid,['%%STEP ' num2str(s) '\n']);
                    fprintf(fid,[num2str(resultCounter+s-1) '\n']);
                end
                if options.addDynamicScalars
                    resultCounter = resultCounter + noStepsUpdated;
                else
                    resultCounter = resultCounter + 1;
                end
                glviewScalarCounter = glviewScalarCounter + 1;
            end
        end

        %% Displacement magnitude (scalars)
        if options.plotMagnitudeDisplacement
            fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
            fprintf(fid,'%%NAME "Magnitude of displacement"\n');
            if options.addDynamicScalars
                noStepsUpdated = noSteps;
            else
                noStepsUpdated = 1;
            end
            for s = 1:noStepsUpdated
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            if options.addDynamicScalars
                resultCounter = resultCounter + noStepsUpdated;
            else
                resultCounter = resultCounter + 1;
            end
            glviewScalarCounter = glviewScalarCounter + 1;
        end

        %% Displacement polar radial direction (scalars)
        if options.plotPolarRadialDisplacement
            fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
            fprintf(fid,'%%NAME "Displacement polar radial direction"\n');
            if options.addDynamicScalars
                noStepsUpdated = noSteps;
            else
                noStepsUpdated = 1;
            end
            for s = 1:noStepsUpdated
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            if options.addDynamicScalars
                resultCounter = resultCounter + noStepsUpdated;
            else
                resultCounter = resultCounter + 1;
            end
            glviewScalarCounter = glviewScalarCounter + 1;
        end

        %% Displacement spherical radial direction (scalars)
        if options.plotSphericalRadialDisplacement
            fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
            fprintf(fid,'%%NAME "Displacement spherical radial direction"\n');
            if options.addDynamicScalars
                noStepsUpdated = noSteps;
            else
                noStepsUpdated = 1;
            end
            for s = 1:noStepsUpdated
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            if options.addDynamicScalars
                resultCounter = resultCounter + noStepsUpdated;
            else
                resultCounter = resultCounter + 1;
            end
            glviewScalarCounter = glviewScalarCounter + 1;
        end

        %% Stress (scalars)
        stressTitles = {'Stress xx','Stress yy','Stress zz','Stress yz','Stress xz','Stress xy'};
        for d = 1:6
            if toggleStressVector(d)
                fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
                fprintf(fid,['%%NAME "' stressTitles{d} '"\n']);
                if options.addDynamicScalars
                    noStepsUpdated = noSteps;
                else
                    noStepsUpdated = 1;
                end
                for s = 1:noStepsUpdated
                    fprintf(fid,['%%STEP ' num2str(s) '\n']);
                    fprintf(fid,[num2str(resultCounter+s-1) '\n']);
                end
                if options.addDynamicScalars
                    resultCounter = resultCounter + noStepsUpdated;
                else
                    resultCounter = resultCounter + 1;
                end
                glviewScalarCounter = glviewScalarCounter + 1;
            end
        end

        %% Polar stress (scalars)
        stressTitles = {'Stress rr','Stress thetatheta','Stress zz','Stress thetaz','Stress rz','Stress rtheta'};
        for d = 1:6
            if togglePolarStressVector(d)
                fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
                fprintf(fid,['%%NAME "' stressTitles{d} '"\n']);
                if options.addDynamicScalars
                    noStepsUpdated = noSteps;
                else
                    noStepsUpdated = 1;
                end
                for s = 1:noStepsUpdated
                    fprintf(fid,['%%STEP ' num2str(s) '\n']);
                    fprintf(fid,[num2str(resultCounter+s-1) '\n']);
                end
                if options.addDynamicScalars
                    resultCounter = resultCounter + noStepsUpdated;
                else
                    resultCounter = resultCounter + 1;
                end
                glviewScalarCounter = glviewScalarCounter + 1;
            end
        end

        %% von Mises Stress (scalars)
        if options.plotVonMisesStress
            fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
            fprintf(fid,'%%NAME "von Mises Stress"\n');
            if options.addDynamicScalars
                noStepsUpdated = noSteps;
            else
                noStepsUpdated = 1;
            end
            for s = 1:noStepsUpdated
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            if options.addDynamicScalars
                resultCounter = resultCounter + noStepsUpdated;
            else
                resultCounter = resultCounter + 1;
            end
            glviewScalarCounter = glviewScalarCounter + 1;
        end
        %% Jacobian (scalars)
        if options.plotOrgJacobi
            fprintf(fid,['\n*GLVIEWSCALAR ' num2str(glviewScalarCounter) '\n']);
            fprintf(fid,'%%NAME "Jacobian"\n');
            if options.addDynamicScalars
                noStepsUpdated = noSteps;
            else
                noStepsUpdated = 1;
            end
            for s = 1:noStepsUpdated
                fprintf(fid,['%%STEP ' num2str(s) '\n']);
                fprintf(fid,[num2str(resultCounter+s-1) '\n']);
            end
            if options.addDynamicScalars
                resultCounter = resultCounter + noStepsUpdated;
            else
                resultCounter = resultCounter + 1;
            end
            glviewScalarCounter = glviewScalarCounter + 1;

        end
    case 'vtk'
        fprintf(fid,'\t\t\t</PointData>\n');
        fprintf(fid,'\t\t</Piece>\n');
        fprintf(fid,'\t</UnstructuredGrid>\n');
        fprintf(fid,'</VTKFile>\n');        
end



fclose(fid);

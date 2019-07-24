function makeVTKfile(data, newOptions)

options = struct('name','untitled',...
                 'celltype','VTK_HEXAHEDRON',...
                 'noSteps',30,...
                 'plotDisplacementVectors',0,...
                 'plotScalarField',0,...
                 'plotFarField',0,...
                 'plotFarFieldError',0,...
                 'plotPolarRadialDisplacement',0,...
                 'plotSphericalRadialDisplacement',0,...
                 'plotStressXX',0,...
                 'plotStressYY',0,...
                 'plotStressZZ',0,...
                 'plotStressYZ',0,...
                 'plotStressXZ',0,...
                 'plotStressXY',0,...
                 'plotdu_xdx',0,...
                 'plotdu_xdy',0,...
                 'plotdu_xdz',0,...
                 'plotdu_ydx',0,...
                 'plotdu_ydy',0,...
                 'plotdu_ydz',0,...
                 'plotdu_zdx',0,...
                 'plotdu_zdy',0,...
                 'plotdu_zdz',0,...
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
                 'plotErrorGrad', 0,...
                 'plotErrorTot', 0,...
                 'plotErrorEnergy', 0,...
                 'plotGradientVectors', 0, ...
                 'plotP_inc', 0, ...
                 'plotAnalytic', 0, ...
                 'plotTotField',0,...
                 'plotTotFieldAbs',0,...
                 'plotSPL',0,...
                 'plotErrorFunc', 0, ...
                 'plotTestFun', 0, ...
                 'plotTestField', 0, ...
                 'timeStepMult', 1, ...
                 'N', NaN, ...
                 'addDynamicScalars', 0);

newOptionFields = fieldnames(newOptions);
for j = 1:numel(newOptionFields)
	options.(newOptionFields{j}) = newOptions.(newOptionFields{j});
end

vtfFileName = options.name;
nodes = data.nodes;
visElements = data.visElements;
if isfield(data,'omega')
    omega = data.omega;
else
    omega = 1;
end
noSteps = options.noSteps;

if ~options.plotTimeOscillation
    noSteps = 1;
end
if length(omega) > 1
    timeStepMult = options.timeStepMult;
    N = options.N;
    T = options.T;
    
    t_arr = linspace(0, T, N+1);
    Nq = N*timeStepMult;
    tq_arr = linspace(0, T*(Nq-1)/Nq, Nq);
    
    %% Displacement (vectors)
    if options.plotDisplacementVectors && timeStepMult > 1
        data.displacement = interp1_new(t_arr, data.displacement, tq_arr);
    end

    %% Gradient (vectors)
    if options.plotGradientVectors && timeStepMult > 1    
        data.gradient = interp1_new(t_arr, data.gradient, tq_arr);
    end

    %% Scalar field
    if options.plotScalarField && timeStepMult > 1
        data.scalarField = interp1_new(t_arr, data.scalarField, tq_arr);
    end

    %% Far field
    if options.plotFarField && timeStepMult > 1
        data.farField = interp1_new(t_arr, data.farField, tq_arr);
    end

    %% Far field Error
    if options.plotFarFieldError && timeStepMult > 1
        data.farFieldError = interp1_new(t_arr, data.farFieldError, tq_arr);
    end
    
    %% Total field
    if options.plotTotField && timeStepMult > 1
        data.totField = interp1_new(t_arr, data.totField, tq_arr);
    end

    %% P_inc
    if options.plotP_inc && timeStepMult > 1
        data.P_inc = interp1_new(t_arr, data.P_inc, tq_arr);
    end


    %% Analytic
    if options.plotAnalytic && timeStepMult > 1
        data.analytic = interp1_new(t_arr, data.analytic, tq_arr);
    end

    %% Error function
    if options.plotErrorFunc && timeStepMult > 1
        data.errorFunc = interp1_new(t_arr, data.errorFunc, tq_arr);
    end

    %% Test function
    if options.plotTestFun && timeStepMult > 1
        data.testFun = interp1_new(t_arr, data.testFun, tq_arr);
    end
    
    %% Test field
    if options.plotTestField && timeStepMult > 1
        data.testField = interp1_new(t_arr, data.testField, tq_arr);
    end
    
    %% Error in scalar field
    if options.plotError && timeStepMult > 1
        data.Error = interp1_new(t_arr, data.Error, tq_arr);
    end
    
    %% Error in gradient field
    if options.plotErrorGrad && timeStepMult > 1
        data.ErrorGrad = interp1_new(t_arr, data.ErrorGrad, tq_arr);
    end
    
    %% Total Error
    if options.plotErrorTot && timeStepMult > 1
        data.ErrorTot = interp1_new(t_arr, data.ErrorTot, tq_arr);
    end

    %% Stress (scalars)
    toggleStressVector = [options.plotStressXX options.plotStressYY options.plotStressZZ ...
                          options.plotStressYZ options.plotStressXZ options.plotStressXY...
                          options.plotPolarStress_rr options.plotPolarStress_thetatheta options.plotPolarStress_zz ...
                            options.plotPolarStress_thetaz options.plotPolarStress_rz options.plotPolarStress_rtheta...
                            options.plotSphericalStress_rr options.plotSphericalStress_thetatheta options.plotSphericalStress_phiphi ...
                                options.plotSphericalStress_thetaphi options.plotSphericalStress_rphi options.plotSphericalStress_rtheta...
                                options.plotVonMisesStress];

    if any(toggleStressVector) && timeStepMult > 1
        data.stress = interp1_new(t_arr, data.stress, tq_arr);
    end
    
    %% Jacobian (scalars)
    toggleJacobianMatrix = [options.plotdu_xdx options.plotdu_xdy options.plotdu_xdz ...
                            options.plotdu_ydx options.plotdu_ydy options.plotdu_ydz ...
                            options.plotdu_zdx options.plotdu_zdy options.plotdu_zdz];

    if any(toggleJacobianMatrix) && timeStepMult > 1
        data.jacobian = interp1_new(t_arr, data.jacobian, tq_arr);
    end
    
else
    Nq = noSteps;
end

[noNodes, d] = size(nodes);
if ~exist(fileparts(vtfFileName),'dir')
    mkdir(fileparts(vtfFileName))
end
if Nq > 1
    delete([vtfFileName '_time_*.vtu'])
end

parfor i = 1:Nq
% for i = 1:Nq
    [noVisElems, noCorners] = size(visElements);
    if Nq == 1
        fid = fopen([vtfFileName '.vtu'],'wt+','b');
    else
        fid = fopen([vtfFileName '_time_' num2str(i) '.vtu'],'wt+','b');
    end
    
    fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid,'\t<UnstructuredGrid>\n');
    fprintf(fid,'\t\t<Piece NumberOfPoints="%d" NumberOfCells="%d">\n',noNodes,noVisElems);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add coordinates
    fprintf(fid,'\t\t\t<Points>\n');
    fprintf(fid,'\t\t\t\t<DataArray type="Float64" NumberOfComponents="%d" format="ascii">\n', d);
    
    for ii = 1:noNodes
        fprintf(fid,['\t\t\t\t\t' repmat('%21.15f ', 1, d) '\n'], nodes(ii,:));
    end
    fprintf(fid,'\t\t\t\t</DataArray>\n');
    fprintf(fid,'\t\t\t</Points>\n');


    fprintf(fid,'\t\t\t<Cells>\n');
    fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    switch options.celltype
        case {'VTK_POLY_LINE', 'VTK_TRIANGLE_STRIP', 'VTK_POLYGON'}
            for eVis = 1:noVisElems
                noCorners = numel(visElements{eVis});
                fprintf(fid,['\t\t\t\t\t' repmat('%d ', 1, noCorners) '\n'], visElements{eVis}-1);
            end            
        otherwise
            for eVis = 1:noVisElems
                fprintf(fid,['\t\t\t\t\t' repmat('%d ', 1, noCorners) '\n'], visElements(eVis,:)-1);
            end
    end
    fprintf(fid,'\t\t\t\t</DataArray>\n');

    fprintf(fid,'\t\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid,'\t\t\t\t\t');
    switch options.celltype
        case {'VTK_POLY_LINE', 'VTK_TRIANGLE_STRIP', 'VTK_POLYGON'}
            counter = 0;
            for eVis = 1:noVisElems
                counter = counter + length(visElements{eVis});
                fprintf(fid,'%d ', counter);
            end
        otherwise
            for eVis = 1:noVisElems
                fprintf(fid,'%d ', noCorners*eVis);
            end
    end
            
    fprintf(fid,'\n');
    fprintf(fid,'\t\t\t\t</DataArray>\n');

    fprintf(fid,'\t\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid,'\t\t\t\t\t');
    switch options.celltype
        case 'VTK_VERTEX'
            fprintf(fid,repmat('1 ', 1, noVisElems));
        case 'VTK_POLY_VERTEX'
            fprintf(fid,repmat('2 ', 1, noVisElems));
        case 'VTK_LINE'
            fprintf(fid,repmat('3 ', 1, noVisElems));
        case 'VTK_POLY_LINE'
            fprintf(fid,repmat('4 ', 1, noVisElems));
        case 'VTK_TRIANGLE'
            fprintf(fid,repmat('5 ', 1, noVisElems));
        case 'VTK_TRIANGLE_STRIP'
            fprintf(fid,repmat('6 ', 1, noVisElems));
        case 'VTK_POLYGON'
            fprintf(fid,repmat('7 ', 1, noVisElems));
        case 'VTK_PIXEL'
            fprintf(fid,repmat('8 ', 1, noVisElems));
        case 'VTK_QUAD'
            fprintf(fid,repmat('9 ', 1, noVisElems));
        case 'VTK_TETRA'
            fprintf(fid,repmat('10 ', 1, noVisElems));
        case 'VTK_VOXEL'
            fprintf(fid,repmat('11 ', 1, noVisElems));
        case 'VTK_HEXAHEDRON'
            fprintf(fid,repmat('12 ', 1, noVisElems));
        case 'VTK_WEDGE'
            fprintf(fid,repmat('13 ', 1, noVisElems));
        case 'VTK_PYRAMID'
            fprintf(fid,repmat('14 ', 1, noVisElems));
        case 'VTK_QUADRATIC_EDGE'
            fprintf(fid,repmat('21 ', 1, noVisElems));
        case 'VTK_QUADRATIC_TRIANGLE'
            fprintf(fid,repmat('22 ', 1, noVisElems));
        case 'VTK_QUADRATIC_QUAD'
            fprintf(fid,repmat('23 ', 1, noVisElems));
        case 'VTK_QUADRATIC_TETRA'
            fprintf(fid,repmat('24 ', 1, noVisElems));
        case 'VTK_QUADRATIC_HEXAHEDRON'
            fprintf(fid,repmat('25 ', 1, noVisElems));
    end
    fprintf(fid,'\n');
    fprintf(fid,'\t\t\t\t</DataArray>\n');
    
    fprintf(fid,'\t\t\t</Cells>\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add results
    fprintf(fid,'\t\t\t<PointData Scalars="scalars">\n');

    %% Displacement (vectors)
    if options.plotDisplacementVectors    
        printField(data.displacement(:,:,i), 'Displacement', fid)
    end

    %% Gradient (vectors)
    if options.plotGradientVectors    
        printField(data.gradient(:,:,i), 'Gradient', fid)
    end

    %% Scalar field
    if options.plotScalarField
        printField(data.scalarField(:,:,i), 'Scalar field', fid)
    end

    %% Far field
    if options.plotFarField
        printField(data.farField(:,:,i), 'Far field', fid)
    end

    %% Far field error
    if options.plotFarFieldError
        printField(data.farFieldError(:,:,i), 'Error far field', fid)
    end

    %% Total field
    if options.plotTotField
        printField(data.totField(:,:,i), 'Total scalar field (real)', fid)
    end

    %% Total field (absolute Value)
    if options.plotTotFieldAbs
        printField(data.totFieldAbs(:,:,i), 'Total scalar field (abs)', fid)
    end

    %% Sound pressure level
    if options.plotSPL
        printField(20*log10(data.totFieldAbs(:,:,i)), 'SPL', fid)
    end

    %% P_inc
    if options.plotP_inc
        printField(data.P_inc(:,:,i), 'P_inc', fid)
    end

    %% Analytic
    if options.plotAnalytic
        printField(data.analytic(:,:,i), 'Analytic', fid)
    end

    %% Error function
    if options.plotErrorFunc
        printField(data.errorFunc(:,:,i), 'plotErrorFunc', fid)
    end

    %% Test function
    if options.plotTestFun
        printField(data.testFun(:,:,i), 'plotTestFun', fid)
    end
    
    %% Test field
    if options.plotTestField
        printField(data.testField(:,:,i), 'plotTestField', fid)
    end
    
    %% Error in scalar field
    if options.plotError
        printField(data.Error(:,:,i), 'Error', fid)
    end
    
    %% Error in gradient field
    if options.plotErrorGrad
        printField(data.ErrorGrad(:,:,i), 'Error Gradient', fid)
    end
    
    %% Total Error
    if options.plotErrorTot
        printField(data.ErrorTot(:,:,i), 'Error Tot', fid)
    end
    
    %% Error in energy
    if options.plotErrorEnergy
        printField(data.ErrorEnergy(:,:,i), 'Error in Energy', fid)
    end

    %% Displacement polar radial direction (scalars)
    if options.plotPolarRadialDisplacement
        if length(omega) > 1
            error('not implemented')
        end
        x = nodes(:,1);
        y = nodes(:,2);
        theta = atan2(y,x);
        displacement = data.displacement;
        u_x = displacement(:,1,i);
        u_y = displacement(:,2,i);
        u_r = u_y.*sin(theta) + u_x.*cos(theta);
        printField(u_r(:,:), 'Polar radial displacement', fid)
    end

    %% Displacement spherical radial direction (scalars)
    if options.plotSphericalRadialDisplacement
        x = nodes(:,1);
        y = nodes(:,2);
        z = nodes(:,3);
        r = sqrt(x.^2+y.^2+z.^2);
        phi = atan2(y,x);
        theta = acos(z./r);
        theta(r < eps) = 0;
        theta = repmat(theta,1,1,Nq);
        phi = repmat(phi,1,1,Nq);
        displacement = data.displacement;
        u_x = displacement(:,1,i);
        u_y = displacement(:,2,i);
        u_z = displacement(:,3,i);
        u_r = u_x.*sin(theta).*cos(phi) + u_y.*sin(theta).*sin(phi) + u_z.*cos(theta);
    %     u_theta = u_x.*cos(theta).*cos(phi) + u_y.*cos(theta).*sin(phi) - u_z.*sin(theta);
    %     u_phi = -u_x.*sin(phi) + u_y.*cos(phi);
        printField(u_r(:,:), 'Spherical radial displacement', fid)
    end

    %% Stress (scalars)
    toggleStressVector = [options.plotStressXX options.plotStressYY options.plotStressZZ ...
                          options.plotStressYZ options.plotStressXZ options.plotStressXY];
    stressTitles = {'Stress xx','Stress yy','Stress zz','Stress yz','Stress xz','Stress xy'};
    
    for l = 1:6
        if toggleStressVector(l)
            stress = data.stress;
            printField(stress(:,l,i), stressTitles{l}, fid)
        end
    end

    %% Derivatives (scalars)
    toggleJacobianMatrix = [options.plotdu_xdx options.plotdu_xdy options.plotdu_xdz ...
                            options.plotdu_ydx options.plotdu_ydy options.plotdu_ydz ...
                            options.plotdu_zdx options.plotdu_zdy options.plotdu_zdz];
    jacobianTitles = {'du_xdx','du_xdy','du_xdz','du_ydx','du_ydy','du_ydz','du_zdx','du_zdy','du_zdz'};
    
    for l = 1:9
        if toggleJacobianMatrix(l)
            jacobian = data.jacobian;
            printField(jacobian(:,l,i), jacobianTitles{l}, fid)
        end
    end
    
    %% Polar stress (scalars)
    togglePolarStressVector = [options.plotPolarStress_rr options.plotPolarStress_thetatheta options.plotPolarStress_zz ...
                                options.plotPolarStress_thetaz options.plotPolarStress_rz options.plotPolarStress_rtheta];
    stressTitles = {'Stress rr','Stress thetatheta','Stress zz','Stress thetaz','Stress rz','Stress rtheta'};
    if any(togglePolarStressVector)
        if length(omega) > 1
            error('not implemented')
        end
        x = nodes(:,1);
        y = nodes(:,2);
        theta = atan2(y,x);
        stress = data.stress;
        polarStress = [cos(theta).^2.*stress(:,1,i) + sin(theta).^2.*stress(:,2,i) + sin(2*theta).*stress(:,6,i), ...
                       sin(theta).^2.*stress(:,1,i) + cos(theta).^2.*stress(:,2,i) - sin(2*theta).*stress(:,6,i), ...
                       stress(:,3,i), ...
                       stress(:,4,i).*cos(theta) - stress(:,5,i).*sin(theta), ...
                       stress(:,4,i).*sin(theta) + stress(:,5,i).*cos(theta), ...
                       -1/2*sin(2*theta).*stress(:,1,i) + 1/2*sin(2*theta).*stress(:,2,i) - cos(2*theta).*stress(:,6,i)];
    end
    for l = 1:6
        if togglePolarStressVector(l)
            printField(polarStress(:,l), stressTitles{l}, fid)
        end
    end

    %% Spherical stress (scalars)
    toggleSphericalStressVector = [options.plotSphericalStress_rr options.plotSphericalStress_thetatheta options.plotSphericalStress_phiphi ...
                                    options.plotSphericalStress_thetaphi options.plotSphericalStress_rphi options.plotSphericalStress_rtheta];
    stressTitles = {'Stress rr','Stress thetatheta','Stress phiphi','Stress thetaphi', 'Stress rphi','Stress rtheta'};
    if any(toggleSphericalStressVector)
        x = nodes(:,1);
        y = nodes(:,2);
        z = nodes(:,3);
        r = sqrt(x.^2+y.^2+z.^2);
        phi = atan2(y,x);
        theta = acos(z./r);
        theta(r < eps) = 0;

        stress = data.stress;
        D = getStressTransformationMatrix(theta,phi,1);
        sphericalStress = zeros(size(stress,1),6);
        for k = 1:6
            for l = 1:6
                D_kl = D(k, l, :);
                sphericalStress(:,k) = sphericalStress(:,k) + D_kl(:).*stress(:,l,i);
            end
        end
    end
    for l = 1:6
        if toggleSphericalStressVector(l)
            printField(sphericalStress(:,l), stressTitles{l}, fid)
        end
    end

    %% von Mises Stress (scalars)
    if options.plotVonMisesStress
        stress = real(data.stress);
        sigma_v = sqrt(((stress(:,1,i)-stress(:,2,i)).^2 + (stress(:,2,i)-stress(:,3,i)).^2 + (stress(:,1,i)-stress(:,3,i)).^2 ...
                                    + 6*(stress(:,4,i).^2 + stress(:,5,i).^2 + stress(:,6,i).^2))/2);
        printField(sigma_v(:,1), 'von Mises Stress', fid)
    end

    %% Jacobian (scalars)
    if options.plotOrgJacobi
        printField(data.orgJacobi(:,:,i), 'Jacobian', fid)
    end

    fprintf(fid,'\t\t\t</PointData>\n');
    fprintf(fid,'\t\t</Piece>\n');
    fprintf(fid,'\t</UnstructuredGrid>\n');
    fprintf(fid,'</VTKFile>');        
    
    fclose(fid);
    if Nq > 1
        disp(['Completed ' num2str(i) ' out of ' num2str(Nq) ' time steps.'])
    end
end

% uncomment the following to get .pvd collector files

% fid = fopen([vtfFileName '.pvd'],'wt+','b');
% 
% fprintf(fid,'<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n');
% fprintf(fid,'\t<Collection>\n');
% for s = 1:noSteps
% %     t = (s-1)/noSteps*2*pi/omega;
% %     fprintf(fid,'\t\t<DataSet timestep="%15.10g" file=''%s\\%s_time_%d.vtu''/>\n',t, cd, vtfFileName, s);
%     fprintf(fid,'\t\t<DataSet timestep="%12.10f" file=''%s\\%s_time_%d.vtu''/>\n',2*(s-1)/noSteps, cd, vtfFileName, s);
% %     fprintf(fid,'\t\t<DataSet timestep="%d" file=''%s\\%s_time_%d.vtu''/>\n',s, cd, vtfFileName, s);
% end
% fprintf(fid,'\t</Collection>\n');
% fprintf(fid,'</VTKFile>');
% 
% fclose(fid);

function printField(field, name, fid)
[m, d] = size(field);

fprintf(fid,'\t\t\t\t<DataArray type="Float64" Name="%s" NumberOfComponents="%d" format="ascii">\n', name, d);
for j = 1:m
    fprintf(fid,['\t\t\t\t\t' repmat(' %21.15g',1,d) '\n'], field(j,:));
end 
fprintf(fid,'\t\t\t\t</DataArray>\n');

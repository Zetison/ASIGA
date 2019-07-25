close all
clear all

memoryCheap = 0;
waringMessageId = 'MATLAB:mir_warning_maybe_uninitialized_temporary';
warning('off',waringMessageId)

startMatlabPool

% matlabpool open 2
% myCluster = parcluster('local')
% parpool(myCluster,2)
% delete(myCluster.Jobs)
Eps = 1e-12;
plot3Dgeometry = 0;
plot2Dgeometry = false;
boundaryMethod = false;
plotTimeOscillation = true;
coreMethod = 'IGA';

E = 50e9; % 50e9–90e9
nu = 0.3; % 0.18–0.3
t1 = 0.001;
t2 = 0.0035; % Thickness of base
rho_s = 2579; % 2400 - 2800 
C = zeros(6,6);
C(1:3,1:3) = E/(1+nu)/(1-2*nu)*...
             [1-nu nu nu;
              nu 1-nu nu; 
              nu nu 1-nu];
C(4:6,4:6) = 0.5*E/(1+nu)*eye(3);
runInParallell = 321;



M = 5;
scaling = 1/10;
degreeElevArr = [1 1 1;
                 2 2 2;
                 3 3 3];
% degreeElevArr = [1 1 1];

stringShift = 40;
for degreeCase = 1:size(degreeElevArr,1)
    fprintf('\n\nCalculating data for case %d\n', degreeCase)
    nurbs = getWineGlassData2();
    Etas = unique(nurbs{2}.knots{2});
    arcLengths = zeros(numel(Etas)-1,1);
    for i = 1:numel(Etas)-1
        arcLengths(i) = arcLength(nurbs{2},Etas(i),Etas(i+1),[0,1],'eta');
    end
    nEta = round(2^(M-1)*arcLengths/scaling);
    iEta = [];
    for i = 1:numel(Etas)-1
        if ismember(i,[2,3,6])
            iEta = [iEta, linspace2(Etas(i),Etas(i+1),round(2^M*arcLengths(i)/scaling))];
        elseif ismember(i,5)
            iEta = [iEta, linspace2(Etas(i),Etas(i+1),round(2^(M+1)*arcLengths(i)/scaling))];
        else
            iEta = [iEta, linspace2(Etas(i),Etas(i+1),nEta(i))];
        end
    end
    nurbs = elevateDegreeInPatches(nurbs,degreeElevArr(degreeCase,:));
    nurbs = insertKnotsInPatches(nurbs,round(2^(M-1)*0.054/scaling),0,round(2^(M-1)*t2/2/scaling));
    nurbs(1) = insertKnotsInPatches(nurbs(1),0,round(2^(M-1)*pi/2*t2/2/scaling),0);
    nurbs(3) = insertKnotsInPatches(nurbs(3),0,round(2^(M-1)*pi/2*t1/2/scaling),0);
    nurbs{2} = insertKnotsInNURBS(nurbs{2},{[] iEta []});
    if plot3Dgeometry
        fluid = nurbs;
        plotMeshAndGeometry
        return
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Collect varables into varCol, varCol_solid and varCol_fluid_i
    varCol.dimension = 3;     
    varCol.C = C;
    varCol = convertNURBS(nurbs, varCol);   

    %% Convert NURBS data
    varCol = generateIGA3DMesh_new(varCol);

    varCol = findDofsToRemove(varCol);    
    fprintf('\nTotal number of elements = %d', varCol.noElems)
% return
    %% Build stiffness matrix
    tic
    fprintf(['\n%-' num2str(stringShift) 's'], 'Building outer fluid matrix ... ')
    options = {'operator','linearElasticity',...
               'fieldDimension', 3,...
               'buildMassMatrix', 1};
    [A_K, A_M] = buildGlobalMatricesVec(varCol, options);
    fprintf('using %12f seconds.', toc)

    %% Modify stiffness matrix and force vector due to glued nodes and homogeneous dirichlet conditions
    dofsToRemove = varCol.dofsToRemove;
    noDofs = varCol.noDofs;
    A_K(dofsToRemove,:) = [];

    A_K(:,dofsToRemove) = [];     

    A_M(dofsToRemove,:) = [];  % zero out the rows and columns of the K matrix

    A_M(:,dofsToRemove) = [];


    %% SOLVE SYSTEM
    tic
    fprintf(['\n%-' num2str(stringShift) 's'], 'Solving eigen value problem ... ')
    noModes = 50;
    [V,D] = eigs(A_K,A_M,noModes+6,'smallestabs','Tolerance',1e-6);
    fprintf('using %12f seconds.', toc)

    eigenFrequencies = sqrt(diag(D(7:end,7:end))/rho_s)/(2*pi); % experimental: 493.88Hz
    p_xi = nurbs{1}.degree(1);
    p_eta = nurbs{1}.degree(2);
    p_zeta = nurbs{1}.degree(3);
    fileNameApp = ['p' num2str(p_xi) 'q' num2str(p_eta) 'r' num2str(p_zeta)];
    fid = fopen(['../plotData/wineGlass/sortedEigenFrequencies' fileNameApp '.txt'],'wt+','b');
    fprintf(fid,'mode\tEigenFrequency\n');
    for i = 1:noModes
        fprintf(fid,'%d\t%1.15f\n', i+6, eigenFrequencies(i));
    end
    fclose(fid);
    continue
    %% POST-PROCESSING
    if degreeCase == size(degreeElevArr,1)
        for i = [6,4,6,8,10,12,13,14,15,17,19,21]
            fprintf(['\n%-' num2str(stringShift) 's'], ['Post processing vibration ' num2str(i) ' ... ' ])
            varCol.omega = sqrt(D(i+6,i+6));
            %% Add solution to removed nodes
            U = zeros(noDofs,1);
            U(setdiff(1:noDofs, dofsToRemove')) = V(:,i+6);  
            U = addSolutionToRemovedNodes_new(U, varCol);

            U = 0.2*U;
            resultsFolderNameParaview = '../../graphics/paraview/wineGlass';
            if ~exist(resultsFolderNameParaview, 'dir')
                mkdir(resultsFolderNameParaview);
            end              
            tic
            vtfFileName = [resultsFolderNameParaview '/mode' num2str(i) '_' fileNameApp];


            extraXiPts = round(10/2^(M-3)); % .. per element
            extraEtaPts = round(10/2^(M-3)); % .. per element
            extraZetaPts = round(2/2^(M-3)); % .. per element
%             extraXiPts = round(5/2^(M-3)); % .. per element
%             extraEtaPts = round(5/2^(M-3)); % .. per element
%             extraZetaPts = round(1/2^(M-3)); % .. per element
%             extraXiPts = 0; % .. per element
%             extraEtaPts = 0; % .. per element
%             extraZetaPts = 0; % .. per element
    %       
            options = struct('name',vtfFileName, 'celltype', 'VTK_HEXAHEDRON', 'plotTimeOscillation', plotTimeOscillation, ...
                                        'plotErrorGrad', 0, 'plotDisplacementVectors',1, 'plotErrorEnergy', 0, ...                   
                            'plotSphericalStress_rr',0, 'plotError', 0,'plotVonMisesStress',1, ...
                             'plotStressXX',1,...
                             'plotStressYY',1,...
                             'plotStressZZ',1,...
                             'plotStressYZ',1,...
                             'plotStressXZ',1,...
                             'plotStressXY',1); 

            varCol.isOuterDomain = false;
            createParaviewFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options, []);
            createVTKmeshFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options)
        end
    end
end

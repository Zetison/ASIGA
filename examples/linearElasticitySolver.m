clear all
close all

addpath(genpath('../export_fig'))
no2Dpoints = 1000;
startMatlabPool
para.plotResultsInParaview = 1;
clear varCol
model = 'bridge-quadratic';
% model = 'ScordelisLoRoof';
prePlot.plot3Dgeometry = 0;
homeDir = expanduser('~');


for M = 2
    switch model
        case 'bridge-quadratic'
%             nurbs = read_g2([homeDir '/OneDrive/SINTEF/DT/bridge-quadratic/' model '.g2']);
            nurbs = read_g2([homeDir '/OneDrive - SINTEF/SINTEF/DT/bridge-quadratic/' model '.g2']);
            E = 3.0e10;
            nu = 0.2;
            rho = 2500;
            % data from https://van.physics.illinois.edu/qa/listing.php?id=64061&t=how-gravitational-force-varies-at-different-locations-on-earth
%             g_45 = 9.806;
%             g_poles = 9.832;
%             g_equator = 9.780;
%             latitude = 59.145905; % at Kjøkøysund bridge
%             latitude = 63.423967; % at Elgseter bridge
%             latitude = 58.16; % at Elgseter bridge
%             g = g_45-1/2*(g_poles-g_equator)*cos(2*latitude*pi/180); ???
            g = 9.81;
%             noInsPts = 0;
%             noInsPts = 1;
            noInsPts = 3;
            newKnots = (2^(M-1)-1)*[1,1,1];
        case 'ScordelisLoRoof' % M = 6 in Cottrell2006iat_Isogeometric_analysis  
            options.degree = 5*ones(1,3);
%             options.degree = [4,5,4];
%             options.degree = [7,7,7];
            options.degree = [2,2,2];
            nurbs = eval(['get' model 'Data(options)']);
            E = 4.32e8;
            nu = 0.0;
            rho = 1;
            t = 0.25;
%             g = 90;
            g = 90/t;
            noInsPts = 100;
            newKnots = (2^(M-1)-1)*[0,1,1]+[3,0,0];
    end
    varCol.BC = 'None';    % not needed
    varCol.media = 'solid';
    varCol.analyticSolutionExist = false;
    varCol.applyLoad = 'None'; % no inhomogeneous Neumann condition
    varCol.boundaryMethod = false;
    nurbs = insertKnotsInNURBS(nurbs,newKnots);
    force = @(v) rho*repmat([0,0,-g],size(v,1),1);

    varCol.extraGP = 0;
    varCol.force = force;
    varCol.C = elasticityMatrix(E,nu);
    varCol.dimension = 3;
    varCol.nurbs = nurbs;
    varCol.operator = 'linearElasticity';
    varCol.fieldDimension = 3;
    varCol.applyBodyLoading = true;
    varCol.buildMassMatrix = false;
    varCol.buildStiffnessMatrix = true;
    varCol.progressBars = 0;
    
    varCol = convertNURBS(varCol);
    varCol = generateIGAmesh(varCol);
    varCol = findDofsToRemove(varCol);
    noDofs = varCol.noDofs;

    stringShift = 40;
    controlPts = varCol.controlPts;
    Eps = 1e-7;
    if prePlot.plot3Dgeometry
        plotNURBS(nurbs,'color',getColor(1),'resolution',10*ones(1,3),'plotControlPolygon',0)
        axis equal
        view(getView())
        camlight
    end
    switch model
        case 'bridge-quadratic'
            fixedXdofs = find(abs(controlPts(:,3)) < Eps).';        
            fixedYdofs = [find(and(abs(controlPts(:,3)-28) < Eps, controlPts(:,1) < -55+Eps)).', ...
                          find(and(abs(controlPts(:,3)-28) < Eps, controlPts(:,1) > 165-Eps)).', ...
                          fixedXdofs];        
            fixedZdofs = fixedYdofs;
            if prePlot.plot3Dgeometry
                plotControlPts2(varCol,{fixedXdofs,  fixedYdofs,  fixedZdofs},...
                                       {'fixedXdofs','fixedYdofs','fixedZdofs'})
            end
            paraExtraScale = 1;
        case 'ScordelisLoRoof'
            L = 50;
            R = 25;
            phi = 40*pi/180;
            fixedZdofs = find(abs(controlPts(:,2) - L/2) < Eps)';
            fixedYdofs = find(abs(controlPts(:,2)) < Eps)';
            fixedXdofs = union(fixedZdofs, find(abs(controlPts(:,1)) < Eps)');
            midSpanNode = find((abs(controlPts(:,3)-(R+t/2)) < Eps).*(abs(controlPts(:,2)) < Eps).*(abs(controlPts(:,1)) < Eps))';
            if prePlot.plot3Dgeometry
                plotControlPts2(varCol,{fixedXdofs,  fixedYdofs,  fixedZdofs, midSpanNode},...
                                       {'fixedXdofs','fixedYdofs','fixedZdofs','midSpanNode'})
            end
            paraExtraScale = 20;
    end

    varCol.progressBars = true;
    d = varCol.dimension;
    fixeddofs = union(union(fixedXdofs*d-2,fixedYdofs*d-1),fixedZdofs*d); 
    dofsToRemove = varCol.dofsToRemove;
    dofsToRemove = union(fixeddofs,dofsToRemove);

    %% Build stiffness matrix
    fprintf(['\n%-' num2str(stringShift) 's'], 'Building system matrices ... ')
    tic
    task.varCol{1} = varCol;
    task.misc.extraGP = [0,0,0];
    task.analyticSolutionExist = false;
    task.splitExteriorFields = false;
    task.misc.BC = 'SHBC';
    task = buildMatrices(task);
    A = task.varCol{1}.A_K;
    F = task.varCol{1}.FF;
    fprintf('using %12f seconds.', toc)

    %% Modify stiffness matrix and force vector due to glued nodes and homogeneous dirichlet conditions
    A(dofsToRemove,:) = [];
    A(:,dofsToRemove) = [];   
    F(dofsToRemove) = [];  
    fprintf(['\n%-' num2str(stringShift) 's'], 'Solving system of equations ... ') 
    tic
    UU = A\F;
    fprintf('using %12f seconds.', toc)
    U = zeros(noDofs,1);
    U(setdiff(1:noDofs, dofsToRemove'),:) = UU;   
    task.varCol{1}.U = U;
    task.varCol = addSolutionToRemovedNodes(task.varCol);
    para.U = task.varCol{1}.U;
    para.name = [homeDir '/results/ASIGA/' model '/' model '_M' num2str(M)];
    if strcmp(model,'ScordelisLoRoof')
        Uz = U(3:d:noDofs);
        abs(min(Uz))
        0.3024
    end
    if para.plotResultsInParaview

        task.varCol{1}.isOuterDomain = false;
        d = task.varCol{1}.dimension;
        isSolid = d == 3;

        para.celltype = 'VTK_HEXAHEDRON';
        para.plotP_inc = false;
        para.plotScalarField = false;
        para.plotTotField = false; 
        para.plotTotFieldAbs = false; 
        para.plotAnalytic = false; 
        para.plotTimeOscillation = false;
        para.computeGrad = false;
        para.plotVonMisesStress = true;
        para.plotStressXX = true;
        para.plotStressYY = true;
        para.plotStressZZ = true;
        para.plotStressYZ = true;
        para.plotStressXZ = true;
        para.plotStressXY = true;
        para.plotError = false; 
        para.plotArtificialBoundary = false;
        para.plotDisplacementVectors = true;

        para.extraXiPts   = round(paraExtraScale/2^(M-1));  % Extra visualization points in the xi-direction per element
        para.extraEtaPts  = round(paraExtraScale/2^(M-1));  % Extra visualization points in the eta-direction per element
        para.extraZetaPts = round(paraExtraScale/2^(M-1));  % Extra visualization points in the zeta-direction per element

        para.i_MS = 1;
        para.plotMesh = true;
        fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting results in Paraview ... ')
        createParaviewFiles(task, 'U', {U}, 'para_options', para);
        fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting mesh files in Paraview ... ')
        tic
        createVTKmeshFiles(task, 'U', {U}, 'para_options', para)
        fprintf('using %12f seconds.\n', toc)
    end
end


function plotControlPts2(varCol,dofs,labels)

controlPts = varCol.controlPts;
noTypes = numel(dofs);
cmap = [1,0,0;
        1,1,0;
        0,0,1;
        1,0,1;
        0,0,0.5;
        0,1,1;
        0,1,0];

hold on
for i = 1:noTypes
    indices = dofs{i};
    p(i) = plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o',...
                'color',cmap(i,:),'MarkerFaceColor',cmap(i,:), 'MarkerEdgeColor', cmap(i,:));
end
legend(p,labels)
% legend show
end
clear all
close all

addpath(genpath('../export_fig'))
no2Dpoints = 1000;
startMatlabPool
para.plotResultsInParaview = 0;
clear varCol
% model = 'bridge-quadratic';
model = 'ScordelisLoRoof';
prePlot.plot3Dgeometry = false;


for M = 3
    switch model
        case 'bridge-quadratic'
            nurbs = read_g2(['../../../SINTEF/DT/bridge-quadratic/' model '.g2']);
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
        case 'ScordelisLoRoof'        
            options.degree = 5*ones(1,3);
            options.degree = [4,5,4];
%             options.degree = [2,2,2];
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
    end
%                 return
    d = varCol.dimension;
    fixeddofs = union(union(fixedXdofs*d-2,fixedYdofs*d-1),fixedZdofs*d); 
    dofsToRemove = varCol.dofsToRemove;
    dofsToRemove = union(fixeddofs,dofsToRemove);

    %% Build stiffness matrix
    fprintf(['\n%-' num2str(stringShift) 's'], 'Building system matrices ... ')
    tic
    varCol = buildMatrices(varCol);
    A = varCol.A_K;
    F = varCol.F;
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
    U = addSolutionToRemovedNodes_new(U, varCol);
    para.U = U;
    celltype = 'VTK_HEXAHEDRON';
    para.name = ['../../../SINTEF/DT/bridge-quadratic/' model 'ASIGA_' num2str(noInsPts)];
    if strcmp(model,'ScordelisLoRoof')
        Uz = U(3:d:noDofs);
        abs(min(Uz))
        0.3024
        para.extraXiPts = 0;
    else
        para.extraXiPts = noInsPts;
    end
    if para.plotResultsInParaview
        para.extraEtaPts = noInsPts;
        para.extraZetaPts = noInsPts;

        varCol.isOuterDomain = false;
        d = varCol.dimension;
        isSolid = d == 3;

        para.celltype = celltype;
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
        fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting results in Paraview ... ')
        createParaviewFiles(varCol, para);
        fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting mesh files in Paraview ... ')
        tic
        createVTKmeshFiles(varCol, para)
        fprintf('using %12f seconds.\n', toc)
    end
end
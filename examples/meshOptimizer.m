function nurbs = meshOptimizer(nurbs,rigidNodes)






J_1 = getJacobian(R,pts,d_p);


% clear all
% close all
% 
% addpath(genpath('../export_fig'))
% no2Dpoints = 1000;
% startMatlabPool
% para.plotResultsInParaview = 1;
% clear varCol
% model = 'cube';
% prePlot.plot3Dgeometry = 1;
% homeDir = expanduser('~');
% 
% 
% 
% for M = 1
%     switch model
%         case 'cube'
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %% Make Jacobian example
%             close all
%             nurbs = getPrismData();
%             nurbs = elevateNURBSdegree(nurbs,1);
%             nurbs2 = nurbs;
%             s = 0.7;
%             nurbs2{1}.coeffs([1,3],end,end-1,end) = [s,s];
%             nurbs2{1}.coeffs([2,3],end-1,end,end) = [s,s];
%             nurbs2{1}.coeffs([1,2],end,end,end-1) = [s,s];
%             
%             %% We want to make an mesh optimizer that recreates nurbs from nurbs2
% %             nurbs = insertKnotsInNURBS(nurbs,[1,2,3]);
% %             nurbs2 = insertKnotsInNURBS(nurbs2,[1,2,3]);
%             if 0
%                 plotNURBS(nurbs,'plotControlPolygon',true,'plotParmDir',0,'plotJacobian',true,'resolution',[100,100,100])
%                 view([77,17])
%                 axis equal
%                 camlight
%                 colorbar
%                 clim([0,1])
%                 
%                 figure
%                 plotNURBS(nurbs2,'plotControlPolygon',true,'plotParmDir',0,'plotJacobian',true,'resolution',[100,100,100])
%                 view([77,17])
%                 axis equal
%                 camlight
%                 colorbar
%                 clim([0,1])
%             end
%             nurbs = nurbs2;
% 
%             E = 3.0e10;
%             nu = 0.2;
%             rho = 2500;
%             noInsPts = 3;
%             newKnots = (2^(M-1)-1)*[1,1,1];
%             paraExtraScale = 20;
%     end
%     varCol.BC = 'None';    % not needed
%     varCol.media = 'solid';
%     varCol.analyticSolutionExist = false;
%     varCol.applyLoad = 'None'; % no inhomogeneous Neumann condition
%     varCol.boundaryMethod = false;
%     nurbs = insertKnotsInNURBS(nurbs,newKnots);
% 
%     varCol.extraGP = 0;
%     varCol.force = @(v,stopForDebug,J) -maxCompression(v,J,stopForDebug);
%     varCol.C = elasticityMatrix(E,nu);
%     varCol.dimension = 3;
%     varCol.nurbs = nurbs;
%     varCol.operator = 'linearElasticity';
%     varCol.fieldDimension = 3;
%     varCol.applyBodyLoading = true;
%     varCol.buildMassMatrix = false;
%     varCol.progressBars = 0;
%     varCol.buildStiffnessMatrix = true;
%     
%     varCol = convertNURBS(varCol);
%     varCol = generateIGAmesh(varCol);
%     varCol = findDofsToRemove(varCol);
%     noDofs = varCol.noDofs;
% 
%     stringShift = 40;
%     controlPts = varCol.controlPts;
%     Eps = 1e-7;
%     if prePlot.plot3Dgeometry
%         plotNURBS(nurbs,'color',getColor(1),'resolution',30*ones(1,3),'plotControlPolygon',1,'plotJacobian',true)
%         axis equal
%         view(getView())
%         camlight
%         colorbar
%         clim([0,1])
%     end
%     switch model
%         case 'cube'
%             fixedXdofs = find(abs(controlPts(:,1)) < Eps).';        
%             fixedYdofs = find(abs(controlPts(:,2)) < Eps).';        
%             fixedZdofs = find(abs(controlPts(:,3)) < Eps).';    
%             rigiddofs = unique([fixedXdofs,fixedYdofs,fixedZdofs]);
%             if prePlot.plot3Dgeometry
%                 plotControlPts2(varCol,{rigiddofs,  rigiddofs,  rigiddofs},...
%                                        {'fixedXdofs','fixedYdofs','fixedZdofs'})
%             end
%     end
% 
%     varCol.progressBars = true;
%     d = varCol.dimension;
%     fixeddofs = union(union(rigiddofs*d-2,rigiddofs*d-1),rigiddofs*d); 
%     dofsToRemove = varCol.dofsToRemove;
%     dofsToRemove = union(fixeddofs,dofsToRemove);
% 
%     %% Build stiffness matrix
%     fprintf(['\n%-' num2str(stringShift) 's'], 'Building system matrices ... ')
%     tic
%     task.varCol{1} = varCol;
%     task.misc.extraGP = [0,0,0];
%     task.analyticSolutionExist = false;
%     task.splitExteriorFields = false;
%     task.misc.BC = 'SHBC';
%     task = buildMatrices(task);
%     A = task.varCol{1}.A_K;
%     F = task.varCol{1}.FF;
%     fprintf('using %12f seconds.', toc)
% 
%     %% Modify stiffness matrix and force vector due to glued nodes and homogeneous dirichlet conditions
%     A(dofsToRemove,:) = [];
%     A(:,dofsToRemove) = [];   
%     F(dofsToRemove) = [];  
%     fprintf(['\n%-' num2str(stringShift) 's'], 'Solving system of equations ... ') 
%     tic
%     UU = A\F;
%     fprintf('using %12f seconds.', toc)
%     U = zeros(noDofs,1);
%     U(setdiff(1:noDofs, dofsToRemove'),:) = UU;   
%     task.varCol{1}.U = U;
%     task.varCol = addSolutionToRemovedNodes(task.varCol);
%     para.U = task.varCol{1}.U;
%     para.name = [homeDir '/results/ASIGA/' model '/' model '_M' num2str(M)];
%     if para.plotResultsInParaview
% 
%         task.varCol{1}.isOuterDomain = false;
%         d = task.varCol{1}.dimension;
%         isSolid = d == 3;
% 
%         para.celltype = 'VTK_HEXAHEDRON';
%         para.plotP_inc = false;
%         para.plotScalarField = false;
%         para.plotTotField = false; 
%         para.plotTotFieldAbs = false; 
%         para.plotAnalytic = false; 
%         para.plotTimeOscillation = false;
%         para.computeGrad = false;
%         para.plotVonMisesStress = true;
%         para.plotStressXX = true;
%         para.plotStressYY = true;
%         para.plotStressZZ = true;
%         para.plotStressYZ = true;
%         para.plotStressXZ = true;
%         para.plotStressXY = true;
%         para.plotError = false; 
%         para.plotArtificialBoundary = false;
%         para.plotDisplacementVectors = true;
% 
%         para.extraXiPts   = round(paraExtraScale/2^(M-1));  % Extra visualization points in the xi-direction per element
%         para.extraEtaPts  = round(paraExtraScale/2^(M-1));  % Extra visualization points in the eta-direction per element
%         para.extraZetaPts = round(paraExtraScale/2^(M-1));  % Extra visualization points in the zeta-direction per element
% 
%         para.i_MS = 1;
%         para.plotMesh = true;
%         fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting results in Paraview ... ')
%         createParaviewFiles(task, 'U', {U}, 'para_options', para);
%         fprintf(['\n%-' num2str(stringShift) 's'], 'Plotting mesh files in Paraview ... ')
%         tic
%         createVTKmeshFiles(task, 'U', {U}, 'para_options', para)
%         fprintf('using %12f seconds.\n', toc)
%     end
% end
% 
% 
% function plotControlPts2(varCol,dofs,labels)
% 
% controlPts = varCol.controlPts;
% noTypes = numel(dofs);
% cmap = [1,0,0;
%         1,1,0;
%         0,0,1;
%         1,0,1;
%         0,0,0.5;
%         0,1,1;
%         0,1,0];
% 
% hold on
% for i = 1:noTypes
%     indices = dofs{i};
%     p(i) = plot3(controlPts(indices,1),controlPts(indices,2),controlPts(indices,3),'o',...
%                 'color',cmap(i,:),'MarkerFaceColor',cmap(i,:), 'MarkerEdgeColor', cmap(i,:));
% end
% legend(p,labels)
% % legend show
% end
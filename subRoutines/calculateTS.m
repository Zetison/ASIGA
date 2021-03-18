function task = calculateTS(task,runTasksInParallel,stringShift)
if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
end

tic
v = getFarFieldPoints(task);

switch task.misc.method
    case {'IE','ABC','IENSG','BA','BEM','PML'}
        if task.ffp.splineBasedNFPcalc
            [zeta0Nodes, noElems, element, element2, index, pIndex] = meshBoundary(task.varCol{1},'Gamma');
            knotVecs = task.varCol{1}.knotVecs;
            elRange = task.varCol{1}.elRange;
            controlPts = task.varCol{1}.controlPts;
            weights = task.varCol{1}.weights;
            U = task.varCol{1}.U;
            noElems = noElems/3;
            nptsPerEl = round(3000/noElems);
            degree = task.varCol{1}.degree(1:2);
            p_h = complex(zeros(nptsPerEl,noElems));
            v = complex(zeros(nptsPerEl,noElems,3));
    %         for e1 = 1:noElems
            parfor e1 = 1:noElems
                e = 3*(e1-1)+1;
                patch = pIndex(e);
                knots = knotVecs{patch};

                Xi_e = zeros(2,2);
                for i = 1:2
                    Xi_e(i,:) = elRange{i}(index(e,i),:);
                end

                sctr = zeta0Nodes(element(e,:));
                pts = controlPts(sctr,:);
                wgts = weights(zeta0Nodes(element2(e,:)),:); % New
                eta = [Xi_e(2,1)+eps,linspace2(Xi_e(2,1),Xi_e(2,2),nptsPerEl-1)].';
                I = findKnotSpans(degree, [0,eta(1)], knots);
                xi = [zeros(size(eta)),eta];
                R = NURBSbasis(I, xi, degree, knots, wgts);
                p_h(:,e1) = R{1}*U(sctr,:);
                v(:,e1,:) = R{1}*pts;
            end
            p_h = p_h(:);
            v = reshape(v,[],3);
            task.ffp.theta = acos(v(:,3)./norm2(v));
            task.ffp.plotFarField = false;
        else
            p_h = calculateScatteredPressure(task, v, 0);
        end
    case 'MFS'
        p_h = calculateScatteredPressureMFS(task, v);
    case 'KDT'
        switch task.varCol{1}.coreMethod
            case 'linear_FEM'
                noElems = task.varCol{1}.noElems;
                element = task.varCol{1}.element;
                tri = NaN(size(element,1),2,3);
                P = task.varCol{1}.controlPts;
                Eps = 1e2*eps;
                for e = 1:noElems
                    sctr = element(e,:);
                    P1 = P(sctr(1),:);
                    P2 = P(sctr(2),:);
                    P3 = P(sctr(3),:);
                    P4 = P(sctr(4),:);
                    tri_e = NaN(1,2,3);
                    if norm(P1-P2) < Eps
                        tri_e(1,1,:) = element(e,[1,4,3]);
                    elseif norm(P1-P3) < Eps
                        tri_e(1,1,:) = element(e,[1,2,4]);
                    elseif norm(P2-P4) < Eps || norm(P3-P4) < Eps
                        tri_e(1,1,:) = element(e,[1,2,3]);
                    else
                        if norm(P2-P3) > norm(P1-P4)
                            tri_e(1,1,:) = element(e,[1,2,4]);
                            tri_e(1,2,:) = element(e,[1,4,3]);
                        else
                            tri_e(1,1,:) = element(e,[1,2,3]);
                            tri_e(1,2,:) = element(e,[2,4,3]);
                        end                                
                    end
                    tri(e,:,:) = tri_e;
                end
                tri = reshape(tri,size(tri,1)*size(tri,2),3);
                tri(any(isnan(tri),2),:) = [];

                %% Find h_max and store results
                task.varCol{1}.h_max = max([norm2(P(tri(:,1),:)-P(tri(:,2),:)); 
                                    norm2(P(tri(:,1),:)-P(tri(:,3),:)); 
                                    norm2(P(tri(:,2),:)-P(tri(:,3),:))]);
                task.varCol{1}.dofs = size(unique(tri,'rows','stable'),1);
                task.varCol{1}.nepw = lambda(1)./task.varCol{1}.h_max;
                task.varCol{1}.noElems = size(tri,1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     trisurf(tri,P(:,1),P(:,2),P(:,3), 'FaceColor', getColor(1))
%                     view(106,26) % sphere and cube
%                     axis off
%                     axis equal
%                     camlight
%                     ax = gca;               % get the current axis
%                     ax.Clipping = 'off';    % turn clipping off
%                     figureFullScreen(gcf)
% %                     
%                     export_fig(['../../graphics/sphericalShell/trianglesParm2_' num2str(task.varCol{1}.noElems)], '-png', '-transparent', '-r300')
%                     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                p_h = kirchApprTri(tri,P,v,task.varCol{1});
            case 'IGA'
                task.varCol{1}.h_max = findMaxElementDiameter(task.varCol{1}.patches);
                task.varCol{1}.nepw = task.varCol{1}.lambda./task.varCol{1}.h_max;
                task.varCol{1}.dofs = task.varCol{1}.noDofs;
%                     p = calculateScatteredPressureBA(task.varCol{1}, Uc{1}, v, 0, plotFarField);
                p_h = calculateScatteredPressureKDT(task, v);
        end
    case 'RT'
        switch scatteringCase
            case 'MS'
                d_vec = task.varCol{1}.d_vec;
                p_h = zeros(size(d_vec,2),1);
%                     for i = 1:size(d_vec,2) %874%
                noIncDir = size(d_vec,2);
                progressBars = task.varCol{1}.progressBars;
                nProgressStepSize = ceil(noIncDir/1000);
                if progressBars
                    ppm = ParforProgMon('Tracing rays: ', noIncDir, nProgressStepSize);
                else
                    ppm = NaN;
                end
%                     for i = 1:noIncDir
                parfor i = 1:noIncDir
                    if progressBars && mod(i,nProgressStepSize) == 0
                        ppm.increment();
                    end
                    varColTemp2 = task.varCol{1};
                    varColTemp2.d_vec = d_vec(:,i);
%                         tic
                    varColTemp2 = createRays(varColTemp2);
%                         fprintf('\nCreating rays in %12f seconds.', toc)
%                         tic
                    varColTemp2 = traceRays(varColTemp2);    
%                         fprintf('\nTracing rays in %12f seconds.', toc)
%                         tic        
                    p_h(i) = calculateScatteredPressureRT(varColTemp2, v(i,:), task.ffp.plotFarField);
%                         fprintf('\nFar field in %12f seconds.', toc)
                end
            otherwise
                task.varCol{1} = createRays(task.varCol{1});
                task.varCol{1} = traceRays(task.varCol{1});            
                p_h = calculateScatteredPressureRT(task.varCol{1}, v, task.ffp.plotFarField);
        end
end
task.p_h = p_h;
if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
end
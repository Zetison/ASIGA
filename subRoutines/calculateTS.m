function task = calculateTS(varCol,task,runTasksInParallel,stringShift)
plotFarField = task.plotFarField;
if ~runTasksInParallel
    fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
end

tic
v = getFarFieldPoints(task.alpha,task.beta,task.r);

switch task.method
    case {'IE','ABC','IENSG','BA','BEM'}
        p_h = calculateScatteredPressure(varCol, v, 0, plotFarField);
    case 'MFS'
        p_h = calculateScatteredPressureMFS(varCol{1}, v, plotFarField);
    case 'KDT'
        switch varCol{1}.coreMethod
            case 'linear_FEM'
                noElems = varCol{1}.noElems;
                element = varCol{1}.element;
                tri = NaN(size(element,1),2,3);
                P = varCol{1}.controlPts;
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
                varCol{1}.h_max = max([norm2(P(tri(:,1),:)-P(tri(:,2),:)); 
                                    norm2(P(tri(:,1),:)-P(tri(:,3),:)); 
                                    norm2(P(tri(:,2),:)-P(tri(:,3),:))]);
                varCol{1}.dofs = size(unique(tri,'rows','stable'),1);
                varCol{1}.nepw = lambda(1)./varCol{1}.h_max;
                varCol{1}.noElems = size(tri,1);
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
%                     export_fig(['../../graphics/sphericalShell/trianglesParm2_' num2str(varCol{1}.noElems)], '-png', '-transparent', '-r300')
%                     
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                p_h = kirchApprTri(tri,P,v,varCol{1});
            case 'IGA'
                varCol{1}.h_max = findMaxElementDiameter(varCol{1}.patches);
                varCol{1}.nepw = varCol{1}.lambda./varCol{1}.h_max;
                varCol{1}.dofs = varCol{1}.noDofs;
%                     p = calculateScatteredPressureBA(varCol{1}, Uc{1}, v, 0, plotFarField);
                p_h = calculateScatteredPressureKDT(varCol{1}, v, plotFarField);
        end
    case 'RT'
        switch scatteringCase
            case 'MS'
                d_vec = varCol{1}.d_vec;
                p_h = zeros(size(d_vec,2),1);
%                     for i = 1:size(d_vec,2) %874%
                noIncDir = size(d_vec,2);
                progressBars = varCol{1}.progressBars;
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
                    varColTemp2 = varCol{1};
                    varColTemp2.d_vec = d_vec(:,i);
%                         tic
                    varColTemp2 = createRays(varColTemp2);
%                         fprintf('\nCreating rays in %12f seconds.', toc)
%                         tic
                    varColTemp2 = traceRays(varColTemp2);    
%                         fprintf('\nTracing rays in %12f seconds.', toc)
%                         tic        
                    p_h(i) = calculateScatteredPressureRT(varColTemp2, v(i,:), plotFarField);
%                         fprintf('\nFar field in %12f seconds.', toc)
                end
            otherwise
                varCol{1} = createRays(varCol{1});
                varCol{1} = traceRays(varCol{1});            
                p_h = calculateScatteredPressureRT(varCol{1}, v, plotFarField);
        end
end
task.results.p = p_h;
task.results.p_Re = real(p_h);
task.results.p_Im = imag(p_h);
task.results.abs_p = abs(p_h);
task.results.TS = 20*log10(abs(p_h/varCol{1}.P_inc));
if varCol{1}.analyticSolutionExist
    if task.plotFarField
        p_ref = varCol{1}.p_0_(v);
%             p_ref = exactKDT(varCol{1}.k,varCol{1}.P_inc,parms.R_o);
    else
        p_ref = varCol{1}.p_(v);
    end
    task.results.p_ref = p_ref;
    task.results.p_Re_ref = real(p_ref);
    task.results.p_Im_ref = imag(p_ref);
    task.results.abs_p_ref = abs(p_ref);
    task.results.TS_ref = 20*log10(abs(p_ref/varCol{1}.P_inc));

    task.results.error_pAbs = 100*abs(abs(p_ref)-abs(p_h))./abs(p_ref);
    task.results.error_p = 100*abs(p_ref-p_h)./abs(p_ref);
end
if ~runTasksInParallel
    fprintf('using %12f seconds.', toc)
end
function task = calculateTS(task,printLog,stringShift)
if printLog
    fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
end

tic
v = getFarFieldPoints(task);

switch task.misc.method
    case {'IE','ABC','IENSG','BA','BEM','PML'}
        if task.ffp.splineBasedNFPcalc
            warning('This implementation is for testing purposes only!')
            varColBdry = meshBoundary(task.varCol{1},'Gamma');

            zeta0Nodes = varColBdry.nodes;
            noElems = varColBdry.noElems;
            element = varColBdry.element;
            element2 = varColBdry.element2;
            index = varColBdry.index;
            pIndex = varColBdry.pIndex;
    
            knotVecs = varColBdry.knotVecs;
            elRange = varColBdry.elRange;
            controlPts = task.varCol{1}.controlPts;
            weights = task.varCol{1}.weights;
            U = task.varCol{1}.U;
            noElems = noElems/3;
            nptsPerEl = round(size(v,1)/noElems);
            degree = task.varCol{1}.degree(1:2);
            p_h = complex(zeros(nptsPerEl,noElems));
            v = zeros(nptsPerEl,noElems,3);
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
            p_h = [p_h(:); U(end)];
            v = reshape(v,[],3);
            task.ffp.theta = [acos(v(:,3)./norm2(v)).',0];
            v = [v; -v(1,:)];
            task.ffp.beta = pi/2-task.ffp.theta;
            task.ffp.plotFarField = false;
        else
            p_h = calculateScatteredPressure(task, v, 0);
        end
    case 'MFS'
        p_h = calculateScatteredPressureMFS(task, v);
    case 'KDT'
        k = task.misc.omega/task.varCol{1}.c_f;
        lambda = 2*pi./k;
        switch task.misc.coreMethod
            case 'linear_FEM'
                task = kirchApprTri(task,v);
                p_h = task.ffp.p_h;
            case 'IGA'
                task.varCol{1}.h_max = findMaxElementDiameter(task.varCol{1}.patches);
                task.varCol{1}.nepw = lambda./task.varCol{1}.h_max;
                task.dofs = task.varCol{1}.noDofs;
                p_h = calculateScatteredPressureKDT(task, v);
        end
    case 'RT'
        warning('RT:limitations','The ray tracing algorithm is purely experimental. Multiple reflection is not implemented, and works only for the S1 or M3 models')
        switch task.misc.scatteringCase
            case 'MS'
                d_vec = task.d_vec;
                p_h = zeros(size(d_vec,2),1);
                noIncDir = size(d_vec,2);
                progressBars = task.misc.progressBars;
                nProgressStepSize = ceil(noIncDir/1000);
                if progressBars
                    ppm = ParforProgMon('Tracing rays: ', noIncDir, nProgressStepSize);
                else
                    ppm = NaN;
                end
%                 for i = 1:noIncDir
                parfor i = 1:noIncDir
                    if progressBars && mod(i,nProgressStepSize) == 0
                        ppm.increment();
                    end
                    task_cp = task;
                    task_cp.d_vec = d_vec(:,i);
                    task_cp = createRays(task_cp);
                    task_cp = traceRays(task_cp);    
                    p_h(i) = calculateScatteredPressureRT(task_cp, v(i,:));
                end
                task.d_vec = d_vec;
            otherwise
                task = createRays(task);
                task = traceRays(task);            
                p_h = calculateScatteredPressureRT(task, v);
        end
end
task.ffp.v = v;
task.ffp.p_h = p_h;
if printLog
    fprintf('using %12f seconds.', toc)
end
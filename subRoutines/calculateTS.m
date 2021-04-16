function task = calculateTS(task,printLog,stringShift)
if printLog
    fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
end

tic
v = getFarFieldPoints(task);

switch task.misc.method
    case {'IE','ABC','IENSG','BA','BEM','PML'}
        if task.ffp.splineBasedNFPcalc
            noElems = task.varCol{1}.noElems;
            element = task.varCol{1}.element;
            element2 = task.varCol{1}.element2;
            index = task.varCol{1}.index;
            pIndex = task.varCol{1}.pIndex;

            d_p = task.varCol{1}.patches{1}.nurbs.d_p;
            knotVecs = task.varCol{1}.knotVecs;
            elRange = task.varCol{1}.elRange;
            controlPts = task.varCol{1}.controlPts;
            weights = task.varCol{1}.weights;
            U = task.varCol{1}.U;
            degree = task.varCol{1}.degree;
            p_h = complex([]);
            v = [];
            paramPts = task.ffp.paramPts;
            for e = 1:noElems
                patch = pIndex(e);
                if isempty(paramPts{patch})
                    continue
                end
                knots = knotVecs{patch};

                Xi_e = zeros(d_p,2);
                for i = 1:d_p
                    Xi_e(i,:) = elRange{i}(index(e,i),:);
                end

                sctr = element(e,:);
                pts = controlPts(sctr,:);
                wgts = weights(element2(e,:),:); % New
                indices = all(and(Xi_e(:,1).' <= paramPts{patch}(:,1:d_p), paramPts{patch}(:,1:d_p) <= Xi_e(:,2).'),2);
                if any(indices)
                    xi = paramPts{patch}(indices,1:d_p);
                    I = findKnotSpans(degree, xi(1,:), knots);
                    R = NURBSbasis(I, xi, degree, knots, wgts, 0);
                    p_h = [p_h; R{1}*U(sctr,:)];
                    v = [v; R{1}*pts];
                end
            end
            task.ffp.theta = acos(v(:,3)./norm2(v)).';
            task.ffp.alpha = atan2(v(:,2),v(:,1));
            task.ffp.beta = pi/2-task.ffp.theta;
            task.ffp.r = norm2(v);
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
                    try
                        ppm = ParforProgMon('Tracing rays: ', noIncDir, nProgressStepSize);
                    catch
                        progressBars = false;
                        ppm = NaN;
                    end
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

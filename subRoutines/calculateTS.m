function task = calculateTS(task,printLog,stringShift)
if printLog
    fprintf(['\n%-' num2str(stringShift) 's'], 'Computing far-field pattern ... ')
end

tic
v = getFarFieldPoints(task);

switch task.misc.method
    case {'IE','ABC','IENSG','BA','BEM','PML'}
        if task.ffp.splineBasedNFPcalc
            U = task.varCol{1}.U;
            p_h = complex([]);
            v = [];
            paramPts = task.ffp.paramPts;
            for patch = 1:task.varCol{1}.noPatches
                if isempty(paramPts{patch})
                    continue
                end
                noCtrlPtsPatch = task.varCol{1}.noCtrlPtsPatch;
                Uindices = (sum(noCtrlPtsPatch(1:patch-1))+1):sum(noCtrlPtsPatch(1:patch));

                xi = paramPts{patch};
                p_h = [p_h; evaluateNURBSvec(task.varCol{1}.nurbs{patch},xi,0,U(Uindices,:))];
                v = [v; evaluateNURBSvec(task.varCol{1}.nurbs{patch},xi)];
            end
            [v,I] = unique(v,'rows');
            p_h = p_h(I,:);
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

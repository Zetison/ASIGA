function task = calculateErrors(task, printLog, stringShift)

omega = task.misc.omega;
if task.err.calculateSurfaceError
    surfaceError = zeros(1,numel(omega));
end
if task.err.calculateSurfEnrgErr
    energyError = zeros(1,numel(omega));
end
if task.err.calculateVolumeError
    L2Error = zeros(1,numel(omega));
    H1Error = zeros(1,numel(omega));
    H1sError = zeros(1,numel(omega));
end
for i_o = 1:numel(omega)
    task.misc.omega = omega(i_o);
    task = getAnalyticSolutions(task);
        
    if task.varCol{1}.boundaryMethod
        if task.err.calculateSurfaceError
            tic 
            if printLog
                fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
            end
            surfaceError = calcSurfErrorBndryMethodVec(task);
            if printLog
                fprintf('using %12f seconds.', toc)   
                fprintf('\nSurface error = %.16g%%', max(surfaceError))
            end
            task.results.surfaceError(i_o) = surfaceError;
        end
        if task.err.calculateSurfEnrgErr
            tic 
            if printLog
                fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface energy error ... ')
            end
            surfEnergyError = calcEnergyErrorBEM(task,1);
            if printLog
                fprintf('using %12f seconds.', toc)   
                fprintf('\nSurface energy error = %.16g%%', max(surfEnergyError))
            end
            task.results.energyError(i_o) = surfEnergyError;
        end
    else
        if task.err.calculateSurfaceError
            tic 
            if printLog
                fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
            end
            if strcmp(task.misc.coreMethod,'SEM')
                surfaceError = calcSurfErrorSEM(task,1);
            else
                surfaceError = calcSurfErrorVec(task,1);
            end
            if printLog
                fprintf('using %12f seconds.', toc) 
                fprintf('\nSurface error = %.16g%%', max(surfaceError))
            end
        end
        if task.err.calculateVolumeError
            tic 
            if printLog
                fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating error ... ')
            end
            if strcmp(task.misc.coreMethod,'SEM')
                [L2Error, H1Error, H1sError, energyError] = calcErrorSEM(task);
            else
                [L2Error, H1Error, H1sError, energyError] = calcErrorVec(task);
            end
            if printLog
                fprintf('using %12f seconds.', toc)
                fprintf('\nVolume L2-error = %.16g%%', max(L2Error))
                fprintf('\nVolume H1-error = %.16g%%', max(H1Error))
                fprintf('\nVolume Energy-error = %.16g%%', max(energyError))
            end
            energyError(i_o) = energyError;
            L2Error(i_o) = L2Error;
            H1Error(i_o) = H1Error;
            H1sError(i_o) = H1sError;
        end
    end
end
task.misc.omega = omega;
if task.err.calculateSurfaceError
    task.results.surfaceError = surfaceError;
end
if task.err.calculateSurfEnrgErr
    task.results.surfEnergyError = surfEnergyError;
end
if task.err.calculateVolumeError
    task.results.L2Error = L2Error;
    task.results.H1Error = H1Error;
    task.results.H1sError = H1sError;
    task.results.energyError = energyError;
end
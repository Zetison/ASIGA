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
                s = sprintf('Calculating surface error (%d/%d) ... ', i_o, numel(omega));
                fprintf(['\n%-' num2str(stringShift) 's'], s)
            end
            surfaceError_i = calcSurfErrorBndryMethodVec(task,i_o);
            if printLog
                fprintf('using %12f seconds.', toc)   
                fprintf('\nSurface error = %.16g%%', max(surfaceError_i))
            end
            surfaceError(i_o) = surfaceError_i;
        end
        if task.err.calculateSurfEnrgErr
            tic 
            if printLog
                s = sprintf('Calculating surface energy error (%d/%d) ... ', i_o, numel(omega));
                fprintf(['\n%-' num2str(stringShift) 's'], s)
            end
            surfEnergyError_i = calcEnergyErrorBEM(task,i_o);
            if printLog
                fprintf('using %12f seconds.', toc)   
                fprintf('\nSurface energy error = %.16g%%', max(surfEnergyError_i))
            end
            energyError(i_o) = surfEnergyError_i;
        end
    else
        if task.err.calculateSurfaceError
            tic 
            if printLog
                s = sprintf('Calculating surface error (%d/%d) ... ', i_o, numel(omega));
                fprintf(['\n%-' num2str(stringShift) 's'], s)
            end
            if strcmp(task.misc.coreMethod,'SEM')
                surfaceError(i_o) = calcSurfErrorSEM(task,1);
            else
                surfaceError(i_o) = calcSurfErrorVec(task,1,i_o);
            end
            if printLog
                fprintf('using %12f seconds.', toc) 
                fprintf('\nSurface error = %.16g%%', max(surfaceError))
            end
        end
        if task.err.calculateVolumeError
            tic 
            if printLog
                s = sprintf('Calculating error (%d/%d) ... ', i_o, numel(omega));
                fprintf(['\n%-' num2str(stringShift) 's'], s)
            end
            if strcmp(task.misc.coreMethod,'SEM')
                [L2Error_i, H1Error_i, H1sError_i, energyError_i] = calcErrorSEM(task);
            else
                [L2Error_i, H1Error_i, H1sError_i, energyError_i] = calcErrorVec(task,i_o);
            end
            if printLog
                fprintf('using %12f seconds.', toc)
                fprintf('\nVolume L2-error = %.16g%%', max(L2Error_i))
                fprintf('\nVolume H1-error = %.16g%%', max(H1Error_i))
                fprintf('\nVolume Energy-error = %.16g%%', max(energyError_i))
            end
            energyError(i_o) = energyError_i;
            L2Error(i_o) = L2Error_i;
            H1Error(i_o) = H1Error_i;
            H1sError(i_o) = H1sError_i;
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
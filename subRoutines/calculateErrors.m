function [L2Error, H1Error, H1sError, energyError, surfaceError] = calculateErrors(task, printLog, stringShift)
L2Error = NaN;
H1Error = NaN;
H1sError = NaN;
energyError = NaN;
surfaceError = NaN;
if task.varCol{1}.boundaryMethod
    if task.err.calculateSurfaceError
        tic 
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorBndryMethodVec(task);
        if printLog
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface error = %.16g', max(surfaceError))
        end
    end
    if task.err.calculateSurfEnrgErr
        tic 
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface energy error ... ')
        end
        energyError = calcEnergyErrorBEM(task,1);
        if printLog
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface energy error = %.16g', max(energyError))
        end
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
            fprintf('\nSurface error = %.16g', max(surfaceError))
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
            fprintf('\nVolume L2-error = %.16g', max(L2Error))
            fprintf('\nVolume H1-error = %.16g', max(H1Error))
            fprintf('\nVolume Energy-error = %.16g', max(energyError))
        end
    end  
end
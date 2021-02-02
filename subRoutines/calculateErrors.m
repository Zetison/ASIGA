function [L2Error, H1Error, H1sError, energyError, surfaceError] = calculateErrors(task, varCol, printLog, stringShift)
L2Error = NaN;
H1Error = NaN;
H1sError = NaN;
energyError = NaN;
surfaceError = NaN;
if varCol{1}.boundaryMethod
    if task.calculateSurfaceError
        tic 
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        surfaceError = calcSurfErrorBndryMethodVec(varCol, task.LpOrder);
        if printLog
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface error = %.16g', surfaceError)
        end
    end
    if task.calculateSurfEnrgErr
        tic 
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface energy error ... ')
        end
        energyError = calcEnergyErrorBEM(varCol{1});
        if printLog
            fprintf('using %12f seconds.', toc)   
            fprintf('\nSurface energy error = %.16g', energyError)
        end
    end
else
    if task.calculateSurfaceError
        tic 
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating surface error ... ')
        end
        if strcmp(varCol{1}.coreMethod,'SEM')
            surfaceError = calcSurfErrorSEM(varCol{1}, task.LpOrder);
        else
            surfaceError = calcSurfErrorVec(varCol{1}, task.LpOrder);
        end
        if printLog
            fprintf('using %12f seconds.', toc) 
            fprintf('\nSurface error = %.16g', surfaceError)
        end
    end
    if task.calculateVolumeError
        tic 
        if printLog
            fprintf(['\n%-' num2str(stringShift) 's'], 'Calculating error ... ')
        end
        if strcmp(varCol{1}.coreMethod,'SEM')
            [L2Error, H1Error, H1sError, energyError] = calcErrorSEM(varCol{1});
        else
            [L2Error, H1Error, H1sError, energyError] = calcErrorVec(varCol);
        end
        if printLog
            fprintf('using %12f seconds.', toc)
            fprintf('\nVolume L2-error = %.16g', L2Error)
            fprintf('\nVolume H1-error = %.16g', H1Error)
            fprintf('\nVolume Energy-error = %.16g', energyError)
        end
    end  
end